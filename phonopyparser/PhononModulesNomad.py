# Copyright 2016-2018 Fawzi Mohamed, Danio Brambila
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# coding=utf-8
#### phonopy parser written by Hagen-Henrik Kowalski and based on the original work of Joerg Mayer on phonopy-FHI-aims

import numpy as np
import past
import math
import os, logging
import json
from fnmatch import fnmatch
from ase.geometry import cell_to_cellpar, crystal_structure_from_cell
from ase.dft.kpoints import special_paths, parse_path_string
try:
    from ase.dft.kpoints import special_points
except ImportError:
    from ase.dft.kpoints import sc_special_points as special_points
from phonopy.units import *
from phonopy.structure.atoms import Atoms
from phonopy.interface.FHIaims import read_aims, write_aims, read_aims_output
from phonopyparser.con import Control
from phonopy import Phonopy
from phonopy.structure.symmetry import Symmetry
from phonopy.file_IO import write_FORCE_CONSTANTS
from phonopy.harmonic.forces import Forces
from phonopy.harmonic.force_constants import get_force_constants
from phonopy.phonon.band_structure import BandStructure
from nomadcore.unit_conversion.unit_conversion import convert_unit_function
from nomadcore.parser_backend import *
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl

from nomad.parsing.legacy import Backend

AimsFrequencyUnitFactors = { 'cm^-1' : VaspToCm, 'THz' : VaspToTHz, 'meV' : 1E3*VaspToEv }
def get_pretty_print(json_object):
    return json.dumps(json_object, sort_keys=True, indent=4, separators=('"', '\n'))

####Depending on the used scipy version the position of atoms at the border of the supercell can be different
####for example if the x coordinate is supposed dto be 0 it can happen that it is at the other end of the supercell
####clean_position checks for such things
def clean_position(scaled_positions):
        scaled_positions = list(scaled_positions)
        for sp in range(len(scaled_positions)):
                for i in range(len(scaled_positions[sp])):
                        if np.float(np.round(scaled_positions[sp][i],7)) >= 1:
                                #print scaled_positions[sp][i]
                                #print 'A'
                                scaled_positions[sp][i] -= 1.0
                                #print scaled_positions[sp][i]
                        elif scaled_positions[sp][i] <= -1e-5:
                                scaled_positions[sp][i] += 1.0
        scaled_positions = np.array(scaled_positions)
        return scaled_positions

###Just for control function to write FORCE_CONSTANTS

def Write_FORCE_CONSTANTS(phonopy_obj, set_of_forces):
        cells = (("initial cell", phonopy_obj.unitcell,"#"), ("supercell", phonopy_obj.supercell,""))
        Nsuper = phonopy_obj.supercell.get_number_of_atoms()
        forces = []
        for (i, disp) in enumerate(phonopy_obj.displacements):
                atom_number = disp[0]
                displacement = disp[1:]
                forces.append(Forces(atom_number, displacement, set_of_forces[i]))
        force_constants = get_force_constants(forces, phonopy_obj.symmetry, phonopy_obj.supercell)
        write_FORCE_CONSTANTS(force_constants, filename='FORCE_CONSTANTS')
        Hessian = get_force_constants(forces, phonopy_obj.symmetry, phonopy_obj.supercell)
####

#### generate_kPath prepares the path genereated by ASE to be used with
#### the function post_process_band
def generate_kPath_ase(cell, symprec):

    eig_val_max = np.real(np.linalg.eigvals(cell)).max()
    eps = eig_val_max*symprec

    lattice = crystal_structure_from_cell(cell, eps)
    paths = parse_path_string(special_paths[lattice])
    points = special_points[lattice]
    k_points = []
    for p in paths:
        k_points.append([points[k] for k in p])
        for index in range(len(p)):
            if p[index] == 'G':
                p[index] = 'Γ'
    parameters = []
    for h, seg in enumerate(k_points):
        for i, path in enumerate(seg):
            parameter = {}
            parameter['npoints'] = 100
            parameter['startname'] = paths[h][i]
            if i == 0 and len(seg) > 2:
                parameter['kstart'] = path
                parameter['kend'] = seg[i+1]
                parameter['endname'] = paths[h][i+1]
                parameters.append(parameter)
            elif i == (len(seg) - 2):
                parameter['kstart'] = path
                parameter['kend'] = seg[i+1]
                parameter['endname'] = paths[h][i+1]
                parameters.append(parameter)
                break
            else:
                parameter['kstart'] = path
                parameter['kend'] = seg[i+1]
                parameter['endname'] = paths[h][i+1]
                parameters.append(parameter)
    return parameters
####

def Collect_Forces_aims(cell_obj, supercell_matrix, displacement, sym, tol = 1e-6):
        symmetry = Symmetry(cell_obj)
        phonopy_obj = Phonopy(cell_obj,
                                supercell_matrix,
                                distance = displacement,
                                symprec = sym)
        supercells = phonopy_obj.get_supercells_with_displacements()
        directories = []
        digits = int( math.ceil( math.log(len(supercells)+1,10) ) ) + 1
        for i in range(len(supercells)):
                directories.append(("phonopy-FHI-aims-displacement-%0" + str(digits) + "d") % (i+1))
        space_group = phonopy_obj.symmetry.get_international_table()
        set_of_forces = []
        Relative_Path = []
        for directory, supercell in zip(directories, supercells):
                aims_out = os.path.join(directory, directory + ".out")
                if not os.path.isfile(aims_out):
                    logging.warn("!!! file not found: %s" % aims_out)
                    cwd = os.getcwd()
                    con_list = os.listdir(cwd)
                    check_var = False
                    for name in con_list:
                        if fnmatch(name, '*.out') == True:
                                aims_out = '%s/%s' % (directory, name)
                                logging.warn(
                                    "Your file seems to have a wrong name proceeding with %s" % aims_out
                                )
                                check_var = True
                                break
                    if check_var == False:
                        raise Exception("No phonon calculations found")
                    os.chdir("../")
                Relative_Path.append(aims_out)
                supercell_calculated = read_aims_output(aims_out)
                if ( (supercell_calculated.get_number_of_atoms() == supercell.get_number_of_atoms()) and
                     (supercell_calculated.get_atomic_numbers() == supercell.get_atomic_numbers()).all() and
                     (abs(supercell_calculated.get_positions()-supercell.get_positions()) < tol).all() and
                     (abs(supercell_calculated.get_cell()-supercell.get_cell()) < tol).all() ):
                    # read_aims_output reads in forces from FHI-aims output as list structure,
                    # but further processing below requires numpy array
                    forces = np.array(supercell_calculated.get_forces())
                    drift_force = forces.sum(axis=0)
                    for force in forces:
                        force -= drift_force / forces.shape[0]
                    set_of_forces.append(forces)
                elif ( (supercell_calculated.get_number_of_atoms() == supercell.get_number_of_atoms()) and
                     (supercell_calculated.get_atomic_numbers() == supercell.get_atomic_numbers()).all() and
                     (abs(clean_position(supercell_calculated.get_scaled_positions())-clean_position(supercell.get_scaled_positions())) < tol).all() and
                     (abs(supercell_calculated.get_cell()-supercell.get_cell()) < tol).all() ):
                     logging.warn("!!! there seems to be a rounding error")
                     forces = np.array(supercell_calculated.get_forces())
                     drift_force = forces.sum(axis=0)
                     for force in forces:
                        force -= drift_force / forces.shape[0]
                     set_of_forces.append(forces)
                else:
                    raise Exception("calculated varies from expected supercell in FHI-aims output %s" % aims_out)

        return set_of_forces, phonopy_obj, Relative_Path


class Get_Properties():
    def __init__(self,
                hessian = None,
                cell = None,
                positions = None,
                symbols = None,
                SC_matrix = None,
                symmetry_thresh = None,
                displacement = None,
                name = None,
                metaInfoEnv = None,
                parser_info = None,
                ):

        #### restoring units
        convert_Phi = convert_unit_function('joules*meter**-2', 'eV*angstrom**-2')
        convert_angstrom = convert_unit_function('meter', 'angstrom')
        hessian = convert_Phi(hessian)
        cell = convert_angstrom(cell)
        positions = convert_angstrom(positions)
        displacement = convert_angstrom(displacement)
        ####

        #### Constructing phonopy_obj
        cell_obj = Atoms(cell = list(cell), symbols= list(symbols), positions= list(positions))
        scaled_positions = cell_obj.get_scaled_positions()
        phonopy_obj = Phonopy(cell_obj, SC_matrix, distance = displacement, symprec = symmetry_thresh)
        phonopy_obj.set_force_constants(hessian)
        ####

        self.phonopy_obj = phonopy_obj

        #### name of the file where the properties are to be stored in
        self.name = name
        self.metaInfoEnv = metaInfoEnv
        self.parser_info = parser_info


        #### choosing mesh
        num_of_atoms = cell_obj.get_number_of_atoms()
        mesh_density = 2*80**3/num_of_atoms
        power_factor = float(1)/float(3)
        mesh_number = np.round(mesh_density**power_factor)
        logging.info('# proceding with a mesh of %d*%d*%d',mesh_number, mesh_number, mesh_number)
        self.mesh = [mesh_number,mesh_number,mesh_number]
        ####

        #### setting parameters
        self.parameters = generate_kPath_ase(cell, symmetry_thresh)
        ####

        #### getting number of atoms
        self.num_of_atoms = num_of_atoms
        self.num_of_atoms_supercell = phonopy_obj.supercell.get_number_of_atoms()

    def post_process_band(self, frequency_unit_factor, parameters = None, is_eigenvectors=False, lookup_labels=False):

        phonopy_obj = self.phonopy_obj
        if parameters == None:
            parameters = self.parameters
        bands = []
        # Distances calculated in phonopy.band_structure.BandStructure object
        # are based on absolute positions of q-points in reciprocal space
        # as calculated by using the cell which is handed over during instantiation.
        # Fooling that object by handing over a "unit cell" diag(1,1,1) instead clashes
        # with calculation of non-analytical terms.
        # Hence generate appropriate distances and special k-points list based on fractional
        # coordinates in reciprocal space (to keep backwards compatibility with previous
        # FHI-aims phonon implementation).
        bands_distances = []
        distance = 0.0
        bands_special_points = [distance]
        bands_labels = []
        label = parameters[0]["startname"]
        for b in parameters:
            kstart = np.array(b["kstart"])
            kend = np.array(b["kend"])
            npoints = b["npoints"]
            dk = (kend-kstart)/(npoints-1)
            bands.append([(kstart + dk*n) for n in range(npoints)])
            dk_length = np.linalg.norm(dk)
            for n in range(npoints):
                bands_distances.append(distance + dk_length*n)
            distance += dk_length * (npoints-1)
            bands_special_points.append(distance)
            label = [b["startname"], b["endname"]]
            if lookup_labels:
                bands_labels.append(BandStructureLabels.get(label.lower(),label))
            else:
                bands_labels.append(label)
        bs_obj = BandStructure(bands,
                               phonopy_obj.dynamical_matrix,
                               is_eigenvectors=is_eigenvectors,
                               factor=frequency_unit_factor)
        freqs = bs_obj.get_frequencies()
        return np.array(freqs), np.array(bands), np.array(bands_labels)

    def get_dos(self, mesh = None):

            phonopy_obj = self.phonopy_obj
            if mesh == None:
                mesh = self.mesh
            phonopy_obj.set_mesh(mesh, is_gamma_center=True)
            q_points = phonopy_obj.get_mesh()[0]
            phonopy_obj.set_qpoints_phonon(q_points, is_eigenvectors = False)
            frequencies = phonopy_obj.get_qpoints_phonon()[0]
            min_freq = min(np.ravel(frequencies))
            max_freq = max(np.ravel(frequencies)) + max(np.ravel(frequencies))*0.05
            phonopy_obj.set_total_DOS(freq_min= min_freq, freq_max = max_freq, tetrahedron_method = True)
            f,dos = phonopy_obj.get_total_DOS()
            return f, dos

    def get_thermodynamical_properties(self, mesh = None, t_max = None, t_min = None, t_step = None):

            phonopy_obj = self.phonopy_obj
            if t_max == None:
                t_max = 1000
            if t_min == None:
                t_min = 0
            if t_step == None:
                t_step = 10
            if mesh == None:
                mesh = self.mesh
            phonopy_obj.set_mesh(mesh, is_gamma_center=True)
            phonopy_obj.set_thermal_properties(t_step = t_step, t_max = t_max, t_min = t_min)
            T, fe, entropy, cv = phonopy_obj.get_thermal_properties()
            kJmolToEv = 1.0 / EvTokJmol
            fe = fe*kJmolToEv
            JmolToEv = kJmolToEv / 1000
            cv = JmolToEv*cv
            return T, fe, entropy, cv

    def prep_bands(self, Emit, parameters = None):

        #name = self.name
        #metaInfoEnv = self.metaInfoEnv
        #parser_info = self.parser_info

        freqs, bands, bands_labels = self.post_process_band(VaspToTHz)

        #### converting THz to eV
        freqs = freqs*THzToEv
        ####

        #### converting eV to Joules
        eVtoJoules = convert_unit_function('eV', 'joules')
        freqs = eVtoJoules(freqs)
        ####

        #### emitting frequencies
        skBand = Emit.openSection("section_k_band")
        Emit.addValue("band_structure_kind", "vibrational")
        for i in range(len(freqs)):
            freq = np.expand_dims(freqs[i], axis = 0)
            skBands = Emit.openSection("section_k_band_segment")
            Emit.addArrayValues("band_energies", freq)
            Emit.addArrayValues("band_k_points", bands[i])
            Emit.addArrayValues("band_segm_labels", bands_labels[i])
            Emit.closeSection("section_k_band_segment", skBands)
        Emit.closeSection("section_k_band", skBand)
        ####

    def prep_density_of_states(self, Emit, mesh = None):

        #name = self.name
        #metaInfoEnv = self.metaInfoEnv
        #parser_info = self.parser_info

        #### Determening DOS
        f, dos = self.get_dos(mesh)
        ####

        #### To match the shape given in metha data another dimension is added to the array (spin degress of fredom is 1)
        dos = np.expand_dims(dos, axis = 0)
        ####

        #### converting THz to eV to Joules
        eVtoJoules = convert_unit_function('eV', 'joules')
        f = f*THzToEv
        f = eVtoJoules(f)
        ####

        #### emitting density of states
        sDos = Emit.openSection("section_dos")
        Emit.addValue("dos_kind", "vibrational")
        Emit.addArrayValues("dos_values", dos)
        Emit.addArrayValues("dos_energies", f)
        Emit.closeSection("section_dos", sDos)
        ####

    def prep_thermodynamical_properties(self, Emit, sSingleConf, mesh = None, t_max = None, t_min = None, t_step = None):

        #name = self.name
        #metaInfoEnv = self.metaInfoEnv
        #parser_info = self.parser_info
        T, fe, entropy, cv = self.get_thermodynamical_properties(mesh = mesh, t_max = t_max, t_min = t_min, t_step = t_step)

        #### deviding free energy by number of atoms to obtain free energy per atom
        fe = fe/self.num_of_atoms
        ####

        # The thermodynamic properties are reported by phonopy for the base
        # system. Since the values in the metainfo are stored per the referenced
        # system, we need to multiple by the size factor between the base system
        # and the supersystem used in the calculations.
        cv = cv*(self.num_of_atoms_supercell/self.num_of_atoms)

        #### converting units
        eVtoJoules = convert_unit_function('eV', 'joules')
        eVperKtoJoules = convert_unit_function('eV*K**-1', 'joules*K**-1')
        fe = eVtoJoules(fe)
        cv = eVperKtoJoules(cv)
        ####

        #### emitting
        frameSeq = Emit.openSection("section_frame_sequence")
        Emit.addArrayValues("frame_sequence_local_frames_ref", np.array([sSingleConf]))
        sTD = Emit.openSection("section_thermodynamical_properties")
        Emit.addArrayValues("thermodynamical_property_temperature", T)
        Emit.addArrayValues("vibrational_free_energy_at_constant_volume", fe)
        Emit.addArrayValues("thermodynamical_property_heat_capacity_C_v", cv)
        sSamplingM = Emit.openSection("section_sampling_method")
        Emit.addValue("sampling_method", "taylor_expansion")
        Emit.addValue("sampling_method_expansion_order", 2)
        Emit.addValue("frame_sequence_to_sampling_ref", sSamplingM)
        Emit.closeSection("section_thermodynamical_properties", sTD)
        Emit.closeSection("section_sampling_method", sSamplingM)
        Emit.closeSection("section_frame_sequence",frameSeq)

    def prep_ref(self, ref_list, Emit):
        sCalc = Emit.openSection("section_calculation_to_calculation_refs")
        Emit.addValue("calculation_to_calculation_kind", "source_calculation")
        for ref in ref_list:
            Emit.addValue("calculation_to_calculation_external_url", ref)
        Emit.closeSection("section_calculation_to_calculation_refs", sCalc)

    def emit_properties(self, emit = ["bands", "dos", "thermodynamical_properties"], parameters = None, mesh = None, t_max = None, t_min = None, t_step = None):

        #### emit has to be either "bands", "dos", or "thermodynamical_properties" default is all of them

        name = self.name
        metaInfoEnv = self.metaInfoEnv
        parser_info = self.parser_info
        Emit = Backend(metaInfoEnv)
        Emit.startedParsingSession(name, parser_info)
        sRun = Emit.openSection("section_run")
        sSingleConf = Emit.openSection("section_single_configuration_calculation")
        for get in emit:
            if get == "bands":
                self.prep_bands(Emit, parameters)
            if get == "dos":
                self.prep_density_of_states(Emit, mesh)
            if get == "thermodynamical_properties":
                self.prep_thermodynamical_properties(Emit, sSingleConf, mesh, t_max, t_min, t_step)
        Emit.closeSection("section_single_configuration_calculation", sSingleConf)
        Emit.closeSection("section_run", sRun)
        Emit.finishedParsingSession("ParseSuccess", None)

    def prem_emit(self,
                Emit,
                sSingleConf,
                emit = ["bands", "dos", "thermodynamical_properties"],
                parameters = None,
                mesh = None,
                t_max = None,
                t_min = None,
                t_step = None):

        for get in emit:
            if get == "bands":
                self.prep_bands(Emit, parameters)
            if get == "dos":
                self.prep_density_of_states(Emit, mesh)
            if get == "thermodynamical_properties":
                self.prep_thermodynamical_properties(Emit, sSingleConf, mesh, t_max, t_min, t_step)


