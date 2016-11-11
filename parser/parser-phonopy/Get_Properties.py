#### phonopy parser based on the original work of Joerg Mayer on phonopy-FHI-aims 

import numpy as np
from PhononModulesNomad import *
from fnmatch import fnmatch
import sys
import math
import os, logging
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
if __name__ == "main":
    import sys
    name = sys.argv[1]
    #### Reading basic properties from JSON
    with open(name) as FORCES:
            data = json.load(FORCES)
    hessian= np.array(data["sections"]["section_single_configuration_calculation-0"]["hessian_matrix"])
    SC_matrix = np.array(data["sections"]["section_system-1"]["SC_matrix"])
    cell = np.array(data["sections"]["section_system-0"]["simulation_cell"])
    symbols = np.array(data["sections"]["section_system-0"]["atom_labels"])
    positions = np.array(data["sections"]["section_system-0"]["atom_positions"])
    displacement = np.array(data["sections"]["section_method-0"]["x_phonopy_displacement"])
    symmetry_thresh = np.array(data["sections"]["section_method-0"]["x_phonopy_symprec"])
    ####


    #### omitting
    get_properties = get_properties(hessian, cell, positions, symbols, SC_matrix, symmetry_thresh, file_name, metaInfoEnv, parser_info)
    get_properties.omit_properties()
