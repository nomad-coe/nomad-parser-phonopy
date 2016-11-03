import numpy as np
from PhononModulesNomad import *
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
from phonopy.file_IO import parse_FORCE_CONSTANTS
from phonopy.harmonic.forces import Forces
from phonopy.harmonic.force_constants import get_force_constants
from phonopy.units import *
from phonopy.structure.atoms import Atoms
from nomadcore.unit_conversion.unit_conversion import convert_unit_function

#### Reading basic properties from JSON
with open('TEST.json') as TEST:
	data = json.load(TEST)
Phi= np.array(data[0]['flatValues']).reshape(data[0]['valuesShape'])
supercell_matrix = np.array(data[1]['flatValues']).reshape(data[1]['valuesShape'])
cell = np.array(data[2]['flatValues']).reshape(data[2]['valuesShape'])
symbols = np.array(data[3]['flatValues']).reshape(data[3]['valuesShape']) #strictly speaking not necessary to reshape
positions = np.array(data[4]['flatValues']).reshape(data[4]['valuesShape'])
displacement = data[5]['value']
sym = data[6]['value']
####

#### Restoring units
convert_Phi = convert_unit_function('joules*meter**-2', 'eV*angstrom**-2')
convert_angstrom = convert_unit_function('meter', 'angstrom')
Phi = convert_Phi(Phi)
cell = convert_angstrom(cell)
positions = convert_angstrom(positions)
displacement = convert_angstrom(displacement)
####

#### Constructing phonopy_obj
cell_obj = Atoms(cell = list(cell), symbols= list(symbols), positions= list(positions))
scaled_positions = cell_obj.get_scaled_positions()
phonopy_obj = Phonopy(cell_obj, supercell_matrix, distance = displacement, symprec = sym)
phonopy_obj.set_force_constants(Phi)

#### for control
print phonopy_obj.symmetry.get_international_table()
####

#### Determening paths in reciprocal space
structure = pm.Structure(list(cell), list(symbols), scaled_positions)
Kpath = HighSymmKpath(structure, symprec=sym, angle_tolerance=5)
kpointpath = dict(Kpath.kpath)
print kpointpath
#from nomadcore.parser_backend import *
parameters = generate_kPath(kpointpath)
post_process_band(phonopy_obj, parameters, VaspToTHz)
num_of_atoms = cell_obj.get_number_of_atoms()
mesh_density = 2*80**3/num_of_atoms
power_factor = float(1)/float(3)
mesh_number = np.round(mesh_density**power_factor)
print '# proceding with a mesh of %d*%d*%d' % (mesh_number, mesh_number, mesh_number)
mesh = [mesh_number,mesh_number,mesh_number]
get_dos(phonopy_obj, mesh)
#get_thermal_properties(phonopy_obj, mesh)
