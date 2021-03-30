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

# phonopy parser based on the original work of Joerg Mayer on phonopy-FHI-aims

import numpy as np
import math
import os
import logging
from fnmatch import fnmatch
from phonopy import Phonopy

from phonopyparser.FHIaims import read_aims_output

class Control:
    def __init__(self, file=None):
        self.phonon = {}
        self.phonon["supercell"] = []
        self.phonon["displacement"] = 0.001
        self.phonon["symmetry_thresh"] = 1E-6
        self.phonon["frequency_unit"] = "cm^-1"
        self.phonon["nac"] = {}
        if file is None:
            self.file = "control.in"
        else:
            self.file = file
        self.read_phonon()

    def read_phonon(self):
        f = open(self.file, 'r')
        try:
            for line in f:
                if not line:
                    break
                fields = line.split()
                nfields = len(fields)
                if (nfields > 2) and (fields[0] == "phonon"):
                    if (fields[1] == "supercell"):
                        if (nfields >= 11):
                            s = map(int, fields[2:11])
                            s = list(s)
                            Smat = np.array(s).reshape(3, 3)
                        elif (nfields >= 5):
                            s = map(int, fields[2:5])
                            s = list(s)
                            Smat = np.diag(s)
                        else:
                            raise Exception("invalid supercell")
                        det_Smat = np.linalg.det(Smat)
                        if (det_Smat == 0):
                            raise Exception("supercell matrix is not invertible")
                        # this consistency check is not strictly necessary, since int function above used to transform the input
                        # already throws an exception when decimal numbers are encountered
                        # keep for consistency (and future changes to that spot?) nevertheless...
                        elif (abs(det_Smat - round(det_Smat)) > 1.0e-6):
                            raise Exception("determinant of supercell differs from integer by more than numerical tolerance of 1.0e-6")
                        self.phonon["supercell"] = s
                    if (fields[1] == "displacement"):
                        self.phonon["displacement"] = float(fields[2])
                    if (fields[1] == "symmetry_thresh"):
                        self.phonon["symmetry_thresh"] = float(fields[2])
                    if (fields[1] == "frequency_unit"):
                        self.phonon["frequency_unit"] = fields[2]
                    if (fields[1] == "nac") and (len(fields) >= 4):
                        if (len(fields) >= 7):
                            delta = tuple(list(map(float, fields[4:7])))
                        else:
                            delta = (0.1, 0.1, 0.1)
                        if (delta[0] == 0.0) and (delta[1] == 0.0) and (delta[2] == 0.0):
                            raise Exception("evaluation of frequencies with non-analytic corrections must be shifted by delta away from Gamma")
                        parameters = {
                            "file": fields[2], "method": fields[3].lower(), "delta": delta}
                        self.phonon["nac"].update(parameters)

        except Exception as e:
            raise(e)

        # supercell is mandatory for all what follows
        if not self.phonon["supercell"]:
            raise Exception("no supercell specified in %s" % self.file)
        f.close()


def clean_position(scaled_positions):
    scaled_positions = list(scaled_positions)
    for sp in range(len(scaled_positions)):
        for i in range(len(scaled_positions[sp])):
            if np.float(np.round(scaled_positions[sp][i], 7)) >= 1:
                scaled_positions[sp][i] -= 1.0
            elif scaled_positions[sp][i] <= -1e-5:
                scaled_positions[sp][i] += 1.0
    scaled_positions = np.array(scaled_positions)

    return scaled_positions


def read_forces_aims(cell_obj, supercell_matrix, displacement, sym, tol=1e-6, logger=None):
    if logger is None:
        logger = logging

    phonopy_obj = Phonopy(cell_obj, supercell_matrix, symprec=sym)
    phonopy_obj.generate_displacements(distance=displacement)
    supercells = phonopy_obj.get_supercells_with_displacements()
    directories = []
    digits = int(math.ceil(math.log(len(supercells) + 1, 10))) + 1
    for i in range(len(supercells)):
        directories.append(("phonopy-FHI-aims-displacement-%0" + str(digits) + "d") % (i + 1))
    set_of_forces = []
    Relative_Path = []
    for directory, supercell in zip(directories, supercells):
        aims_out = os.path.join(directory, directory + ".out")
        if not os.path.isfile(aims_out):
            logger.warn("!!! file not found: %s" % aims_out)
            cwd = os.getcwd()
            con_list = os.listdir(cwd)
            check_var = False
            for name in con_list:
                if fnmatch(name, '*.out'):
                    aims_out = '%s/%s' % (directory, name)
                    logger.warn(
                        "Your file seems to have a wrong name proceeding with %s" % aims_out
                    )
                    check_var = True
                    break
            if not check_var:
                logger.error("No phonon calculations found")
                return set_of_forces, phonopy_obj, Relative_Path
            os.chdir("../")
        Relative_Path.append(aims_out)
        supercell_calculated = read_aims_output(aims_out)
        if (
            (supercell_calculated.get_number_of_atoms() == supercell.get_number_of_atoms()) and
            (supercell_calculated.get_atomic_numbers() == supercell.get_atomic_numbers()).all() and
            (abs(supercell_calculated.get_positions() - supercell.get_positions()) < tol).all() and
            (abs(supercell_calculated.get_cell() - supercell.get_cell()) < tol).all()):
            # read_aims_output reads in forces from FHI-aims output as list structure,
            # but further processing below requires numpy array
            forces = np.array(supercell_calculated.get_forces())
            drift_force = forces.sum(axis=0)
            for force in forces:
                force -= drift_force / forces.shape[0]
            set_of_forces.append(forces)
        elif (
            (supercell_calculated.get_number_of_atoms() == supercell.get_number_of_atoms()) and
            (supercell_calculated.get_atomic_numbers() == supercell.get_atomic_numbers()).all() and
            (abs(clean_position(supercell_calculated.get_scaled_positions()) - clean_position(supercell.get_scaled_positions())) < tol).all() and
            (abs(supercell_calculated.get_cell() - supercell.get_cell()) < tol).all()):
            logger.warn("!!! there seems to be a rounding error")
            forces = np.array(supercell_calculated.get_forces())
            drift_force = forces.sum(axis=0)
            for force in forces:
                force -= drift_force / forces.shape[0]
            set_of_forces.append(forces)
        else:
            logger.error("calculated varies from expected supercell in FHI-aims output")

    return set_of_forces, phonopy_obj, Relative_Path
