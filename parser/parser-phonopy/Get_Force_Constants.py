import numpy as np
from PhononModulesNomadpy import *
from fnmatch import fnmatch
import sys
import math
import os
import pymatgen as pm
from pymatgen.symmetry.bandstructure import HighSymmKpath
from phonopy.interface.FHIaims import read_aims, write_aims, read_aims_output
from con import Control
from phonopy import Phonopy
from phonopy.structure.symmetry import Symmetry
from phonopy.file_IO import write_FORCE_CONSTANTS
from phonopy.harmonic.forces import Forces
from phonopy.harmonic.force_constants import get_force_constants
from phonopy.units import *
from nomadcore.unit_conversion.unit_conversion import convert_unit_function

#### determening properties of the undisplaced cell
cell_obj = read_aims("geometry.in")
cell = cell_obj.get_cell()
positions = cell_obj.get_positions()
symbols = np.array(cell_obj.get_chemical_symbols())
control = Control()
if (len(control.phonon["supercell"]) == 3):
                supercell_matrix = np.diag(control.phonon["supercell"])
elif (len(control.phonon["supercell"]) == 9):
                supercell_matrix = np.array(control.phonon["supercell"]).reshape(3,3)
displacement = control.phonon["displacement"]
sym = control.phonon["symmetry_thresh"]
####

#### constructing FORCE_CONSTANTS
set_of_forces, phonopy_obj = Collect_Forces_aims(cell_obj, supercell_matrix, displacement, sym)
phonopy_obj.set_forces(set_of_forces)
phonopy_obj.produce_force_constants()
FC2 = phonopy_obj.get_force_constants()
####

#### obtaining information about supercell
super_c=phonopy_obj.supercell
s_cell = super_c.get_cell()
super_pos = super_c.get_positions()
super_sym = np.array(super_c.get_chemical_symbols())
####

#### Converting properties to Si unitis
converter_FC2 = convert_unit_function('eV*angstrom**-2', 'joules*meter**-2')
convert_angstrom = convert_unit_function('angstrom', 'meter')
FC2 = converter_FC2(FC2)
cell = convert_angstrom(cell)
positions = convert_angstrom(positions)
displacement = convert_angstrom(displacement)
####

#### Writing JSON
from nomadcore.parser_backend import *
TEST = open('TEST.json','w')
path = '../../../../nomad-meta-info/meta_info/nomad_meta_info/public.nomadmetainfo.json'
metaInfoPath = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), path))
metaInfoEnv, warns = loadJsonFile(filePath=metaInfoPath,
                                  dependencyLoader=None,
                                  extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                  uri=None)
Parse = JsonParseEventsWriterBackend(metaInfoEnv)
Parse.startedParsingSession(...)
#sRun = Parse.openSection("section_run")
#sMethod = Parse.openSection("section_method")
#sBaseSystem = Parse.openSection("section_system")
Parse.openSection('section_system')
Parse.addArrayValues('atom_labels', symbols)
Parse.addArrayValues('atom_positions', positions)
Parse.addArrayValue("original_system_ref", cell)
Parse.addArrayValues('SC_Matrix', supercell_matrix)
#Parse.closeSection("section_system", sBaseSystem)
#...
#sSuperCellSystem = Parse.openSection("section_system")
Parse.addArrayValues('super_cell_atom_labels', super_sym)
Parse.addArrayValues('atom_positions', super_pos)
Parse.addArrayValue('original_system_ref', s_cell)
#Parse.closeSection("section_system", sSuperCellSystem)
Parse.openSection('')
#Parse.addArray(None,np.shape(FC2))
#sHessian = 
Parse.addArrayValues('Hessian', FC2)
#Parse.addArray(None,np.shape(supercell_matrix))
Parse.addArrayValues('SC_Matrix', supercell_matrix)
Parse.addArrayValues('Primitive_cell', cell)
Parse.addArrayValues('Primitive_atoms', symbols)
Parse.addArrayValues('postions', positions)
Parse.addValue('displacement', displacement)
Parse.addValue('symprec', sym)
####

