import numpy as np
from phonopy.units import VaspToTHz as AimsToTHz, VaspToCm as AimsToCm, VaspToEv as AimsToEv, THzToCm, THzToEv
from phonopy.interface.FHIaims import read_aims, write_aims, read_aims_output
from phonopy.harmonic.dynamical_matrix import get_equivalent_smallest_vectors, DynamicalMatrix, DynamicalMatrixNAC
from phonopy.structure.cells import get_reduced_bases
import numpy as np
from phonopy import Phonopy
from phonopy.structure.cells import Primitive



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
                            Smat = np.array(s).reshape(3,3)
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
                        elif (abs(det_Smat-round(det_Smat)) > 1.0e-6):
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
                        parameters = { "file" : fields[2],
                                       "method" : fields[3].lower(),
                                       "delta" : delta }
                        self.phonon["nac"].update(parameters)

        except Exception:
            #print (line,)
            #print ("|-> line triggered exception: ") + str(Exception)
            raise
        # supercell is mandatory for all what follows
        if not self.phonon["supercell"]:
            raise Exception("no supercell specified in %s" % self.file)
        f.close()
        
