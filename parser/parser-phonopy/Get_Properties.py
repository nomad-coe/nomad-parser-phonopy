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


    #### Determening paths in reciprocal space
    get_properties = get_properties(hessian, cell, positions, symbols, SC_matrix, symmetry_thresh)
    freqs, bands, bands_labels = get_properties.post_process_band(VaspToTHz)
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
    sRun = Parse.openSection("section_run")
    skBand = Parse.openSection("section_k_band")
    for i in range(len(freqs)):
        freq = np.expand_dims(freqs[i], axis = 0)
        skBands = Parse.openSection("section_k_band_segment")
        Parse.addArrayValues("band_energies", freq)
        Parse.addArrayValues("band_k_points", bands[i])
        Parse.addArrayValues("band_segm_labels", bands_labels[i])
        Parse.close("section_k_band_segment", skBands)
    Parse.close("section_k_band", skBand)
    ####

    #### Determening DOS
    f, dos = get_properties.get_dos(phonopy_obj, mesh)
    ####

    #### To match the shape given in metha data another dimension is added to the array (spin degress of fredom is 1)
    dos = np.expand_dims(dos, axis = 0)
    ####

    #### converting eV to Joules
    f = eVtoJoules(f)
    ####

    #### Parsing dos
    sDos = Parse.openSection("section_dos")
    Parse.addArray("dos_values", dos)
    Parse.addArray("dos_energies", f)
    Parse.close("section_dos", sDos)
    ####

    #### Determening Thermal properties
    T, fe, entropy, cv = get_properties.get_thermal_properties()
    ####
    sHarmonic = Parse.openSection("thermodynamical_properties_calculation_method")
    #### Parsing 
    sTD = Parse.openSection("section_thermodynamical_properties")
    Parse.addArrayValues("thermodynamical_property_temperature", T)
    Parse.addArrayValues("vibrational_free_energy_at_constant_volume", fe)
    Parse.addArrayValues("thermodynamical_property_heat_capacity_C_v", cv)
    ####
    Parse.close("thermodynamical_properties_calculation_method", sHarmonic)
    Parse.close("section_run", sRun)
    Parse.finishedParsingSession("ParseSuccess", None)
