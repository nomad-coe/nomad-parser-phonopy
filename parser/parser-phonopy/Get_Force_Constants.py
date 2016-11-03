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
TEST.write('[')
Parse = JsonParseEventsWriterBackend('TEST',TEST)
Parse.startedParsingSession(...)
sRun = Parse.openSection("section_run")
sMethod = Parse.openSection("section_method")
sBaseSystem = Parse.openSection("section_system")
# output base geometry
Parse.addArrayValues('simulation_cell', cell)
Parse.addArrayValues('atom_labels', symbols)
Parse.addArrayValues('atom_positions', positions)
Parse.closeSection("section_system", sBaseSystem)
sSuperCellSystem = Parse.openSection("section_system")
Parse.addValue("original_system_ref", sBaseSystem)
...
Parse.closeSection("section_system", sSuperCellSystem)

#Parse.addArray(None,np.shape(FC2))
Parse.addArrayValues('Hessian', FC2)
#Parse.addArray(None,np.shape(supercell_matrix))
Parse.addArrayValues('SC_Matrix', supercell_matrix)
Parse.addArrayValues('Primitive_cell', cell)
Parse.addArrayValues('Primitive_atoms', symbols)
Parse.addArrayValues('postions', positions)
Parse.addValue('displacement', displacement)
Parse.addValue('symprec', sym)
TEST.write(']')
####

