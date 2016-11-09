import numpy as np
import past
import math
import os
import json
from fnmatch import fnmatch
from phonopy.units import *
from phonopy.interface.FHIaims import read_aims, write_aims, read_aims_output
from con import Control
from phonopy import Phonopy
from phonopy.structure.symmetry import Symmetry
from phonopy.file_IO import write_FORCE_CONSTANTS
from phonopy.harmonic.forces import Forces
from phonopy.harmonic.force_constants import get_force_constants
from phonopy.phonon.band_structure import BandStructure
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


#### generate_kPath prepares the path genereated by pymatgen to be used with 
#### the function post_process_band
def generate_kPath(kpointpath):
        parameters = []
        if float(np.shape(kpointpath['path'])>1):
                for i, paths in enumerate(kpointpath['path']):
                        for j, path in enumerate(paths):
                                parameter = {}
                                parameter['npoints'] = 100
                                parameter['startname'] = str(path).split()[0]
                                if j < len(paths)-1:
                                        parameter['endname'] = str(paths[j+1]).split()[0]
                                        parameter['kend'] = list(kpointpath['kpoints'][paths[j+1]])
                                        #print kpointpath['kpoints'][path[j+1]]
                                else:
                                        parameter['endname'] = str(kpointpath['path'][0][0]).split()[0]
                                        parameter['kend'] = list(kpointpath['kpoints'][kpointpath['path'][0][0]])
                                parameter['kstart'] = list(kpointpath['kpoints'][path])
                                parameters.append(parameter)

        return parameters

def Collect_Forces_aims(cell_obj, supercell_matrix, displacement, sym, tol = 1e-6):
        symmetry = Symmetry(cell_obj)
        #control = Control()
        phonopy_obj = Phonopy(cell_obj, 
                                supercell_matrix, 
                                distance = displacement, 
                                symprec = sym)
        supercells = phonopy_obj.get_supercells_with_displacements()
        directories = []
        digits = int( math.ceil( math.log(len(supercells)+1,10) ) ) + 1
        for i in range(len(supercells)):
                directories.append( ("phonopy-FHI-aims-displacement-%0" + str(digits) + "d") % (i+1))
        space_group = phonopy_obj.symmetry.get_international_table()
        print (space_group)
        set_of_forces = []
        for directory, supercell in zip(directories, supercells):
                aims_out = os.path.join(directory, directory + ".out")
                if not os.path.isfile(aims_out):
                    print ("!!! file not found: %s" % aims_out)
                    os.chdir(directory)
                    cwd = os.getcwd()
                    con_list = os.listdir(cwd)
                    check_var = False
                    for name in con_list:
                        if fnmatch(name, '*.out') == True:
                                aims_out = '%s/%s' % (directory, name)
                                print ("!!! WARNING your file seems to have a wrong name proceeding with %s" % aims_out)
                                check_var = True
                                break
                    if check_var == False:
                        print ("!!! No phonon calculations found")
                        sys.exit(1)
                    os.chdir("../")
                supercell_calculated = read_aims_output(aims_out)
                if ( (supercell_calculated.get_number_of_atoms() == supercell.get_number_of_atoms()) and
                     (supercell_calculated.get_atomic_numbers() == supercell.get_atomic_numbers()).all() and
                     (abs(supercell_calculated.get_positions()-supercell.get_positions()) < tol).all() and
                     (abs(supercell_calculated.get_cell()-supercell.get_cell()) < tol).all() ):
                    # read_aims_output reads in forces from FHI-aims output as list structure, 
                    # but further processing below requires numpy array
                    forces = np.array(supercell_calculated.get_forces())
                    drift_force = forces.sum(axis=0)
                    #print ("#   | correcting drift : %11.5f %11.5f %11.5f" % tuple(drift_force))
                    for force in forces:
                        force -= drift_force / forces.shape[0]
                    set_of_forces.append(forces)
                elif ( (supercell_calculated.get_number_of_atoms() == supercell.get_number_of_atoms()) and
                     (supercell_calculated.get_atomic_numbers() == supercell.get_atomic_numbers()).all() and
                     (abs(clean_position(supercell_calculated.get_scaled_positions())-clean_position(supercell.get_scaled_positions())) < tol).all() and
                     (abs(supercell_calculated.get_cell()-supercell.get_cell()) < tol).all() ):
                     print ("!!! there seems to be a rounding error")
                     forces = np.array(supercell_calculated.get_forces())
                     drift_force = forces.sum(axis=0)
                     #print "#   | correcting drift : %11.5f %11.5f %11.5f" % tuple(drift_force)
                     for force in forces:
                        force -= drift_force / forces.shape[0]
                     set_of_forces.append(forces)
                else:
                    print ("!!! calculated varies from expected supercell in FHI-aims output %s" % aims_out)
                    sys.exit(2)
        return set_of_forces, phonopy_obj



def post_process_band(phonopy_obj, parameters, frequency_unit_factor, is_eigenvectors=False, write_yaml=False, do_matplotlib=False, lookup_labels=False):
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
        # one long list to simplify output
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
    return freqs, np.array(bands), np.array(bands_labels)

def get_dos(phonopy_obj, mesh):
        phonopy_obj.set_mesh(mesh, is_gamma_center=True)
        q_points = phonopy_obj.get_mesh()[0]
        phonopy_obj.set_qpoints_phonon(q_points, is_eigenvectors = False)
        frequencies = phonopy_obj.get_qpoints_phonon()[0]
        min_freq = min(np.ravel(frequencies))
        max_freq = max(np.ravel(frequencies)) + max(np.ravel(frequencies))*0.05
        phonopy_obj.set_total_DOS(freq_min= min_freq, freq_max = max_freq, tetrahedron_method = True)
        f,dos = phonopy_obj.get_total_DOS()
        return f, dos

def get_thermal_properties(phonopy_obj, mesh):
        print ('#### NOT THE REAL UNITS')
        phonopy_obj.set_mesh(mesh, is_gamma_center=True)
        phonopy_obj.set_thermal_properties()
        T, fe, entropy, cv = phonopy_obj.get_thermal_properties()
        kJmolToEv = 1.0 / EvTokJmol
        JmolToEv = kJmolToEv / 1000
        return T, fe, entropy, cv
