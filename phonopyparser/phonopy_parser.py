#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import os
import numpy as np
import logging
import phonopy
from phonopy.units import THzToEv
from phonopy.structure.atoms import PhonopyAtoms

from phonopyparser.phonopy_properties import PhononProperties

import nomad.config
from nomad.units import ureg
from nomad.parsing.file_parser import TextParser, Quantity
from nomad.datamodel.metainfo.common_dft import Run, System, SystemToSystemRefs, Method,\
    SingleConfigurationCalculation, KBand, KBandSegment, Dos, FrameSequence,\
    ThermodynamicalProperties, SamplingMethod, Workflow, Phonon, CalculationToCalculationRefs

from nomad.parsing.parser import FairdiParser
from .metainfo import m_env


def read_aims(filename):
    '''Method to read FHI-aims geometry files in phonopy context.'''
    cell = []
    positions = []
    fractional = []
    symbols = []
    magmoms = []
    with open(filename) as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.split()
            if line[0] == 'lattice_vector':
                cell.append([float(x) for x in line[1:4]])
            elif line[0].startswith('atom'):
                fractional.append(line[0] == 'atom_frac')
                positions.append([float(x) for x in line[1:4]])
                symbols.append(line[4])
            elif line[0] == 'initial_moment':
                magmoms.append(float(line[1]))

    for n, pos in enumerate(positions):
        if fractional[n]:
            positions[n] = [sum([pos[j] * cell[j][i] for j in range(3)]) for i in range(3)]
    if len(magmoms) == len(positions):
        return PhonopyAtoms(cell=cell, symbols=symbols, positions=positions, magmoms=magmoms)
    else:
        return PhonopyAtoms(cell=cell, symbols=symbols, positions=positions)


class Atoms_with_forces(PhonopyAtoms):
    ''' Hack to phonopy.atoms to maintain ASE compatibility also for forces.'''

    def get_forces(self):
        return self.forces


def read_aims_output(filename):
    ''' Read FHI-aims output
        returns geometry with forces from last self-consistency iteration'''
    cell = []
    symbols = []
    positions = []
    forces = []
    N = 0

    with open(filename) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if 'Number of atoms' in line:
                N = int(line.split()[5])
            elif '| Unit cell:' in line:
                cell = [[float(x) for x in f.readline().split()[1:4]] for _ in range(3)]
            elif 'Atomic structure:' in line or 'Updated atomic structure:' in line:
                positions = []
                symbols = []
                symbol_index = 3 if 'Atomic' in line else 4
                position_index = 4 if 'Atomic' in line else 1
                while len(positions) != N:
                    line = f.readline()
                    if 'Species' in line or 'atom ' in line:
                        line = line.split()
                        positions.append([float(x) for x in line[position_index:position_index + 3]])
                        symbols.append(line[symbol_index])
            elif 'Total atomic forces' in line:
                forces = [[float(x) for x in f.readline().split()[2:5]] for _ in range(N)]

    atoms = Atoms_with_forces(cell=cell, symbols=symbols, positions=positions)
    atoms.forces = forces

    return atoms


def read_forces_aims(reference_supercells, tolerance=1E-6, logger=None):
    '''
    Collect the pre calculated forces for each of the supercells
    '''
    def get_aims_output_file(directory):
        files = [f for f in os.listdir(directory) if f.endswith('.out')]
        output = None
        for f in files:
            try:
                output = read_aims_output(os.path.join(directory, f))
                break
            except Exception:
                pass
        return output

    def is_equal(reference, calculated):
        if len(reference) != len(calculated):
            logger.warn('Inconsistent number of atoms.')
            return False
        if (reference.get_atomic_numbers() != calculated.get_atomic_numbers()).any():
            logger.warn('Inconsistent species.')
            return False
        if (abs(reference.get_cell() - calculated.get_cell()) > tolerance).any():
            logger.warn('Inconsistent cell.')
            return False
        ref_pos = reference.get_scaled_positions()
        cal_pos = calculated.get_scaled_positions()
        # wrap to bounding cell
        ref_pos %= 1.0
        cal_pos %= 1.0
        if (abs(ref_pos - cal_pos) > tolerance).any():
            logger.warn('Inconsistent positions.')
            return False
        return True

    reference_paths, forces_sets = [], []

    n_pad = int(np.ceil(np.log10(len(reference_supercells) + 1))) + 1
    for n, reference_supercell in enumerate(reference_supercells):
        directory = 'phonopy-FHI-aims-displacement-%s' % (str(n + 1).zfill(n_pad))
        filename = os.path.join(directory, '%s.out' % directory)
        if os.path.isfile(filename):
            calculated_supercell = read_aims_output(filename)
        else:
            # try reading out files
            calculated_supercell = get_aims_output_file(directory)

        # compare if calculated cell really corresponds to supercell
        if not is_equal(reference_supercell, calculated_supercell):
            logger.error('Supercells do  not match')

        forces = np.array(calculated_supercell.get_forces())
        drift_force = forces.sum(axis=0)
        for force in forces:
            force -= drift_force / forces.shape[0]
        forces_sets.append(forces)
        reference_paths.append(filename)
    return forces_sets, reference_paths


class ControlParser(TextParser):
    def __init__(self):
        super().__init__()

    def init_quantities(self):
        def str_to_nac(val_in):
            val = val_in.strip().split()
            nac = dict(file=val[0], method=val[1].lower())
            if len(val) > 2:
                nac['delta'] = [float(v) for v in val[3:6]]
            return nac

        def str_to_supercell(val_in):
            val = [int(v) for v in val_in.strip().split()]
            if len(val) == 3:
                return np.diag(val)
            else:
                return np.reshape(val, (3, 3))

        self._quantities = [
            Quantity('displacement', r'\n *phonon displacement\s*([\d\.]+)', dtype=float),
            Quantity('symmetry_thresh', r'\n *phonon symmetry_thresh\s*([\d\.]+)', dtype=float),
            Quantity('frequency_unit', r'\n *phonon frequency_unit\s*(\S+)'),
            Quantity('supercell', r'\n *phonon supercell\s*(.+)', str_operation=str_to_supercell),
            Quantity('nac', r'\n *phonon nac\s*(.+)', str_operation=str_to_nac)]


class PhonopyParser(FairdiParser):
    def __init__(self, **kwargs):
        super().__init__(
            name='parsers/phonopy', code_name='Phonopy', code_homepage='https://phonopy.github.io/phonopy/',
            mainfile_name_re=(r'(.*/phonopy-FHI-aims-displacement-0*1/control.in$)|(.*/phon.+yaml)')
        )
        self._kwargs = kwargs
        self.control_parser = ControlParser()

    @property
    def mainfile(self):
        return self._filepath

    @mainfile.setter
    def mainfile(self, val):
        self._phonopy_obj = None
        self.references = []
        self._filepath = os.path.abspath(val)

    @property
    def calculator(self):
        if 'control.in' in self.mainfile:
            return 'fhi-aims'
        elif self.mainfile.endswith('.yaml'):
            return 'vasp'

    @property
    def phonopy_obj(self):
        if self._phonopy_obj is None:
            if self.calculator == 'fhi-aims':
                self._build_phonopy_object_fhi_aims()
            elif self.calculator == 'vasp':
                self._build_phonopy_object_vasp()
        return self._phonopy_obj

    def _build_phonopy_object_vasp(self):
        cwd = os.getcwd()
        os.chdir(os.path.dirname(self.mainfile))

        try:
            phonopy_obj = phonopy.load(self.mainfile)
        finally:
            os.chdir(cwd)

        self._phonopy_obj = phonopy_obj

    def _build_phonopy_object_fhi_aims(self):
        cwd = os.getcwd()
        os.chdir(os.path.dirname(os.path.dirname(self.mainfile)))

        try:
            cell_obj = read_aims('geometry.in')
            self.control_parser.mainfile = 'control.in'
            supercell_matrix = self.control_parser.get('supercell')
            displacement = self.control_parser.get('displacement', 0.001)
            sym = self.control_parser.get('symmetry_thresh', 1e-6)
            try:
                phonopy_obj = phonopy.Phonopy(cell_obj, supercell_matrix, symprec=sym)
                phonopy_obj.generate_displacements(distance=displacement)
                supercells = phonopy_obj.get_supercells_with_displacements()
                set_of_forces, relative_paths = read_forces_aims(supercells, logger=self.logger)
            except Exception:
                self.logger.error("Error generating phonopy object.")
                set_of_forces = []
                phonopy_obj = None
                relative_paths = []

            prep_path = self.mainfile.split("phonopy-FHI-aims-displacement-")
            # Try to resolve references as paths relative to the upload root.
            try:
                for path in relative_paths:
                    abs_path = "%s%s" % (prep_path[0], path)
                    rel_path = abs_path.split(nomad.config.fs.staging + "/")[1].split("/", 3)[3]
                    self.references.append(rel_path)
            except Exception:
                self.logger.warn("Could not resolve path to a referenced calculation within the upload.")

        finally:
            os.chdir(cwd)

        if set_of_forces:
            phonopy_obj.set_forces(set_of_forces)
            phonopy_obj.produce_force_constants()

        self._phonopy_obj = phonopy_obj

    def parse_bandstructure(self):
        freqs, bands, bands_labels = self.properties.get_bandstructure()
        if freqs is None:
            return

        # convert THz to eV
        freqs = freqs * THzToEv

        # convert eV to J
        freqs = (freqs * ureg.eV).to('joules').magnitude

        sec_scc = self.archive.section_run[0].section_single_configuration_calculation[0]

        sec_k_band = sec_scc.m_create(KBand)
        sec_k_band.band_structure_kind = 'vibrational'

        for i in range(len(freqs)):
            freq = np.expand_dims(freqs[i], axis=0)
            sec_k_band_segment = sec_k_band.m_create(KBandSegment)
            sec_k_band_segment.band_energies = freq
            sec_k_band_segment.band_k_points = bands[i]
            sec_k_band_segment.band_segm_labels = bands_labels[i]

    def parse_dos(self):
        f, dos = self.properties.get_dos()

        # To match the shape given in meta data another dimension is added to the
        # array (spin degress of fredom is 1)
        dos = np.expand_dims(dos, axis=0)

        # convert THz to eV to Joules
        f = f * THzToEv
        f = (f * ureg.eV).to('joules').magnitude

        sec_scc = self.archive.section_run[0].section_single_configuration_calculation[0]
        sec_dos = sec_scc.m_create(Dos)
        sec_dos.dos_kind = 'vibrational'
        sec_dos.dos_values = dos
        sec_dos.dos_energies = f

    def parse_thermodynamical_properties(self):
        T, fe, _, cv = self.properties.get_thermodynamical_properties()

        n_atoms = len(self.phonopy_obj.unitcell)
        n_atoms_supercell = len(self.phonopy_obj.supercell)

        fe = fe / n_atoms

        # The thermodynamic properties are reported by phonopy for the base
        # system. Since the values in the metainfo are stored per the referenced
        # system, we need to multiple by the size factor between the base system
        # and the supersystem used in the calculations.
        cv = cv * (n_atoms_supercell / n_atoms)

        # convert to SI units
        fe = (fe * ureg.eV).to('joules').magnitude

        cv = (cv * ureg.eV / ureg.K).to('joules/K').magnitude

        sec_run = self.archive.section_run[0]
        sec_scc = sec_run.section_single_configuration_calculation

        sec_frame_sequence = sec_run.m_create(FrameSequence)
        sec_frame_sequence.frame_sequence_local_frames_ref = sec_scc

        sec_thermo_prop = sec_frame_sequence.m_create(ThermodynamicalProperties)
        sec_thermo_prop.thermodynamical_property_temperature = T
        sec_thermo_prop.vibrational_free_energy_at_constant_volume = fe
        sec_thermo_prop.thermodynamical_property_heat_capacity_C_v = cv

        sec_sampling_method = sec_run.m_create(SamplingMethod)
        sec_sampling_method.sampling_method = 'taylor_expansion'
        sec_sampling_method.sampling_method_expansion_order = 2

        sec_frame_sequence.frame_sequence_to_sampling_ref = sec_sampling_method

    def parse_ref(self):
        sec_scc = self.archive.section_run[0].section_single_configuration_calculation[0]
        for ref in self.references:
            sec_calc_refs = sec_scc.m_create(CalculationToCalculationRefs)
            sec_calc_refs.calculation_to_calculation_kind = 'source_calculation'
            sec_calc_refs.calculation_to_calculation_external_url = ref

    def parse(self, filepath, archive, logger, **kwargs):
        self.mainfile = os.path.abspath(filepath)
        self.archive = archive
        self.logger = logger if logger is not None else logging
        self._kwargs.update(kwargs)

        self._metainfo_env = m_env

        sec_run = self.archive.m_create(Run)
        sec_run.program_name = 'Phonopy'
        sec_run.program_version = phonopy.__version__

        phonopy_obj = self.phonopy_obj

        pbc = np.array((1, 1, 1), bool)

        unit_cell = self.phonopy_obj.unitcell.get_cell()
        unit_pos = self.phonopy_obj.unitcell.get_positions()
        unit_sym = np.array(self.phonopy_obj.unitcell.get_chemical_symbols())

        super_cell = self.phonopy_obj.supercell.get_cell()
        super_pos = self.phonopy_obj.supercell.get_positions()
        super_sym = np.array(self.phonopy_obj.supercell.get_chemical_symbols())

        unit_cell = (unit_cell * ureg.angstrom).to('meter').magnitude
        unit_pos = (unit_pos * ureg.angstrom).to('meter').magnitude

        super_cell = (super_cell * ureg.angstrom).to('meter').magnitude
        super_pos = (super_pos * ureg.angstrom).to('meter').magnitude

        try:
            displacement = np.linalg.norm(phonopy_obj.displacements[0][1:])
            displacement = (displacement * ureg.angstrom).to('meter').magnitude
        except Exception:
            displacement = None

        supercell_matrix = phonopy_obj.supercell_matrix
        sym_tol = phonopy_obj.symmetry.tolerance

        sec_system_unit = sec_run.m_create(System)
        sec_system_unit.configuration_periodic_dimensions = pbc
        sec_system_unit.atom_labels = unit_sym
        sec_system_unit.atom_positions = unit_pos
        sec_system_unit.simulation_cell = unit_cell

        sec_system = sec_run.m_create(System)
        sec_system_to_system_refs = sec_system.m_create(SystemToSystemRefs)
        sec_system_to_system_refs.system_to_system_kind = 'subsystem'
        sec_system_to_system_refs.system_to_system_ref = sec_system_unit
        sec_system.configuration_periodic_dimensions = pbc
        sec_system.atom_labels = super_sym
        sec_system.atom_positions = super_pos
        sec_system.simulation_cell = super_cell
        sec_system.SC_matrix = supercell_matrix
        sec_system.x_phonopy_original_system_ref = sec_system_unit

        sec_method = sec_run.m_create(Method)
        # TODO I put this so as to have a recognizable section method, but metainfo
        # should be expanded to include phonon related method parameters
        sec_method.electronic_structure_method = 'DFT'
        sec_method.x_phonopy_symprec = sym_tol
        if displacement is not None:
            sec_method.x_phonopy_displacement = displacement

        try:
            force_constants = phonopy_obj.get_force_constants()
            force_constants = (force_constants * ureg.eV / ureg.angstrom ** 2).to('J/(m**2)').magnitude
        except Exception:
            self.logger.error('Error producing force constants.')
            return

        sec_scc = sec_run.m_create(SingleConfigurationCalculation)
        sec_scc.single_configuration_calculation_to_system_ref = sec_system
        sec_scc.single_configuration_to_calculation_method_ref = sec_method
        sec_scc.hessian_matrix = force_constants

        # get bandstructure configuration file
        maindir = os.path.dirname(self.mainfile)
        files = [f for f in os.listdir(maindir) if f.endswith('.conf')]
        self._kwargs.update({'band_conf': os.path.join(maindir, files[0]) if files else None})
        self.properties = PhononProperties(self.phonopy_obj, self.logger, **self._kwargs)

        self.parse_bandstructure()
        self.parse_dos()
        self.parse_thermodynamical_properties()
        self.parse_ref()

        sec_workflow = self.archive.m_create(Workflow)
        sec_workflow.workflow_type = 'phonon'
        sec_phonon = sec_workflow.m_create(Phonon)
        sec_phonon.force_calculator = self.calculator
        vol = np.dot(unit_cell[0], np.cross(unit_cell[1], unit_cell[2]))
        sec_phonon.mesh_density = np.prod(self.properties.mesh) / vol
        n_imaginary = np.count_nonzero(self.properties.frequencies < 0)
        sec_phonon.n_imaginary_frequencies = n_imaginary
        if phonopy_obj.nac_params:
            sec_phonon.with_non_analytic_correction = True
