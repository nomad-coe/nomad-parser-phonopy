#### phonopy parser based on the original work of Joerg Mayer on the phonopy-FHI-aims code

import numpy as np
from PhononModulesNomad import *
from fnmatch import fnmatch
import sys
import math
import os
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

parser_info = {"name": "parser_phonopy", "version": "0.1"}

path = "../../../../nomad-meta-info/meta_info/nomad_meta_info/phonopy.nomadmetainfo.json"
metaInfoPath = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), path))
metaInfoEnv, warns = loadJsonFile(filePath=metaInfoPath,
                                  dependencyLoader=None,
                                  extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                  uri=None)

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
####

#### Determening paths in reciprocal space
parameters = generate_kPath_ase(cell)
freqs, bands, bands_labels = post_process_band(phonopy_obj, parameters, VaspToTHz)
####

#### converting THz to eV
freqs = freqs*THzToEv
####

#### converting eV to Joules
eVtoJoules = convert_unit_function('eV', 'joules')
freqs = eVtoJoules(freqs)
####

#### Parsing frequencies
Parse = JsonParseEventsWriterBackend(metaInfoEnv)
Parse.startedParsingSession(name, parser_info)
skBand = Parse.openSection("section_k_band")
for i in range(len(freqs)):
    freq = np.expand_dims(freqs[i], axis = 0)
    skSegment = Parse.openSection("section_k_band_segment")
    Parse.addArrayValues(freq, skSegment)
    Parse.close("section_k_band_segment", skSegmment)
    skPoints = Parse.openSection("section_k_band_segment")
    Parse.addArrayValues(bands, skPoints)
    Parse.close("section_k_band_segment", skPoints)
    skLabels = Parse.openSection("section_k_band_segment")
    Parse.addArrayValues(bands_labels[i], skLabels)
    Parse.close("section_k_band_segment", skLabels)
Parse.close("section_k_band", skBand)
####

#### Determening DOS
num_of_atoms = cell_obj.get_number_of_atoms()
mesh_density = 2*80**3/num_of_atoms
power_factor = float(1)/float(3)
mesh_number = np.round(mesh_density**power_factor)
print ('# proceding with a mesh of %d*%d*%d' % (mesh_number, mesh_number, mesh_number))
mesh = [mesh_number,mesh_number,mesh_number]
f, dos = get_dos(phonopy_obj, mesh)
####

#### To match the shape given in metha data another dimension is added to the array (spin degress of fredom is 1)
dos = np.expand_dims(dos, axis = 0)
####

#### converting eV to Joules
f = eVtoJoules(f)
####

#### Parsing dos
sDos = Parse.openSection("section_dos")
Parse.addArray(dos, sDos)
Parse.addArray(f, sDos)
Parse.close("section_dos", sDos)
####

#### Determening Thermal properties
T, fe, entropy, cv = get_thermal_properties(phonopy_obj, mesh)
####
