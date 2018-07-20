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

#### phonopy parser written by Hagen-Henrik Kowalski and based on the original work of Joerg Mayer on phonopy-FHI-aims

import numpy as np
from PhononModulesNomad import *
import setup_paths
from fnmatch import fnmatch
import sys
import math
import os
import argparse
from phonopy.interface.FHIaims import read_aims, write_aims, read_aims_output
from con import Control
from phonopy import Phonopy, __version__
from phonopy.structure.symmetry import Symmetry
from phonopy.file_IO import write_FORCE_CONSTANTS
from phonopy.harmonic.forces import Forces
from phonopy.harmonic.force_constants import get_force_constants
from phonopy.units import *
from nomadcore.unit_conversion.unit_conversion import convert_unit_function
from nomadcore.parser_backend import *
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl

phonopy_version = __version__
parser_info = {"name": "parser_phonopy", "version": "1.0"}

path = "../../../../nomad-meta-info/meta_info/nomad_meta_info/phonopy.nomadmetainfo.json"
metaInfoPath = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), path))
metaInfoEnv, warns = loadJsonFile(filePath=metaInfoPath,
                                  dependencyLoader=None,
                                  extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                  uri=None)


def parse(name):
    pbc = np.array((1, 1, 1), bool)
    Parse = JsonParseEventsWriterBackend(metaInfoEnv)
    Parse.startedParsingSession(name, parser_info)
    sRun = Parse.openSection("section_run")
    Parse.addValue("program_name", "Phonopy")
    Parse.addValue("program_version", phonopy_version)
    Basesystem = Parse.openSection("section_system")
    Parse.addArrayValues("configuration_periodic_dimensions", pbc)
    Parse.addArrayValues("atom_labels", symbols)
    Parse.addArrayValues("atom_positions", positions)
    Parse.addArrayValues("simulation_cell", cell)
    Parse.closeSection("section_system", Basesystem)
    Supercellsystem = Parse.openSection("section_system")
    sysrefs = Parse.openSection("section_system_to_system_refs")
    Parse.addValue("system_to_system_kind", "subsystem")
    Parse.addValue("system_to_system_ref", Basesystem)
    Parse.closeSection("section_system_to_system_refs", sysrefs)
    Parse.addArrayValues("configuration_periodic_dimensions", pbc)
    Parse.addArrayValues("atom_labels", super_sym)
    Parse.addArrayValues("atom_positions", super_pos)
    Parse.addArrayValues("simulation_cell", s_cell)
    Parse.addArrayValues("SC_matrix", supercell_matrix)
    Parse.addValue("x_phonopy_original_system_ref", Basesystem)
    Parse.closeSection("section_system", Supercellsystem)
    method = Parse.openSection("section_method")
    Parse.addValue("x_phonopy_symprec", sym)
    Parse.addValue("x_phonopy_displacement", displacement)
    Parse.closeSection("section_method", method)
    results = Parse.openSection("section_single_configuration_calculation")
    Parse.addValue("single_configuration_calculation_to_system_ref", Supercellsystem)
    Parse.addValue("single_configuration_to_calculation_method_ref", method)
    Parse.addArrayValues("hessian_matrix", FC2)
    GP = Get_Properties(FC2, cell, positions, symbols, supercell_matrix, sym, displacement)
    GP.prem_emit(Parse, results)
    GP.prep_ref(Whole_Path, Parse)
    Parse.closeSection("section_single_configuration_calculation", results)
    Parse.closeSection("section_run", sRun)
    Parse.finishedParsingSession("ParseSuccess", None)

#### determening properties of the undisplaced cell
if __name__ == '__main__':
    import sys

    parser = argparse.ArgumentParser(description='Parses a phonopy calculation.')
    parser.add_argument('mainFileUri',
                        help='The uri of the main file associated with this calculation.')
    parser.add_argument('mainFilePath', default = None,
                        help='The path to the main file associated with this calculation.')
    parser.add_argument('--kind', default = 'FHI-aims', choices = ["FHI-aims"],
                        help='The kind of phonopy calculation performed')

    args = parser.parse_args()
    if args.mainFilePath:
        mainDir = os.path.dirname(os.path.dirname(os.path.abspath(args.mainFilePath)))
        os.chdir(mainDir)
    name = args.mainFileUri
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
    set_of_forces, phonopy_obj, Relative_Path = Collect_Forces_aims(cell_obj, supercell_matrix, displacement, sym)
    Prep_Path = name.split("phonopy-FHI-aims-displacement-")
    Whole_Path = []
    for Path in Relative_Path:
        Whole_Path.append("%s%s" % (Prep_Path[0], Path))
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
    s_cell = convert_angstrom(s_cell)
    super_pos = convert_angstrom(super_pos)
    positions = convert_angstrom(positions)
    displacement = convert_angstrom(displacement)
    
    #### parsing
    parse(name)

