import numpy as np
from phonopy.units import VaspToTHz as AimsToTHz, VaspToCm as AimsToCm, VaspToEv as AimsToEv, THzToCm, THzToEv
from phonopy.interface.FHIaims import read_aims, write_aims, read_aims_output
from phonopy.harmonic.dynamical_matrix import get_equivalent_smallest_vectors, DynamicalMatrix, DynamicalMatrixNAC
from phonopy.structure.cells import get_reduced_bases
import numpy as np
from phonopy import Phonopy
from phonopy.structure.cells import Primitive


class needed_methods():
	def get_smallest_vectors(super_c, primitive, symprec=1e-5):

	    p2s_map = primitive.get_primitive_to_supercell_map()
	    size_super = super_c.get_number_of_atoms()
	    size_prim = primitive.get_number_of_atoms()
	    shortest_vectors = np.zeros((size_super, size_prim, 27, 3), dtype='double')
	    multiplicity = np.zeros((size_super, size_prim), dtype='intc')
	
	    for i in range(size_super): # run in supercell
	        for j, s_j in enumerate(p2s_map): # run in primitive
	            vectors = get_equivalent_smallest_vectors(i,
	                                                      s_j,
	                                                      super_c,
	                                                      primitive.get_cell(),
	                                                      symprec)
	            multiplicity[i][j] = len(vectors)
	            for k, elem in enumerate(vectors):
	                shortest_vectors[i][j][k] = elem
	    return shortest_vectors, multiplicity

	def back_to_reality(shortest_vec,super_c,primitive):

	        dim = 3
	        size_super = super_c.get_number_of_atoms()
	        size_prim = primitive.get_number_of_atoms()
	        lattice_p = primitive.get_cell()
	        lattice = np.array(lattice_p)
	        real_vec = np.zeros((size_super,size_prim,27,dim),dtype='double')
	        for i in range(size_super):
	                for j in range(size_prim):
	                        for k in range(27):
	                                real_vec[i][j][k] = np.dot(lattice,shortest_vec[i][j][k])
	
	        return real_vec
		


class Control:
    def __init__(self, file=None):
        self.phonon = {}
        self.phonon["supercell"] = []
        self.phonon["displacement"] = 0.001
        self.phonon["symmetry_thresh"] = 1E-6
        self.phonon["frequency_unit"] = "cm^-1"
        self.phonon["nac"] = {}
        self.phonon["hessian"] = []
        self.phonon["band"] = []
        self.phonon["dos"] = {}
        self.phonon["free_energy"] = {}
        self.phonon["animation"] = []
        self.phonon["modulations"] = {}
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
                    if (fields[1] == "hessian"):
                        self.phonon["hessian"] = fields[2:]
                    if (fields[1] == "band") and (len(fields) >= 11):
                        parameters = { "kstart" : list(map(float, fields[2:5])),
                                       "kend" : list(map(float, fields[5:8])),
                                       "npoints" : int(fields[8]),
                                       "startname" : fields[9],
                                       "endname" : fields[10] }
                        self.phonon["band"].append(parameters)
                    if (fields[1] == "dos") and (len(fields) >= 7):
                        parameters = { "fstart" : float(fields[2]),
                                       "fend" : float(fields[3]),
                                       "fpoints" : int(fields[4]),
                                       "broad" : float(fields[5]),
                                       "qdensity" : list(map(int, fields[6:])) }
                        self.phonon["dos"].update(parameters)
                    if (fields[1] == "free_energy") and (len(fields) >= 6):
                        parameters = { "Tstart" : float(fields[2]),
                                       "Tend" : float(fields[3]),
                                       "Tpoints" : int(fields[4]),
                                       "qdensity" : list(map(int, fields[5:])) }
                        self.phonon["free_energy"].update(parameters)
                    if (fields[1] == "animation") and (len(fields) >= 12):
                        parameters = { "q" : list(map(float, fields[2:5])),
                                       "band" : int(fields[5]),
                                       "amp" : float(fields[6]),
                                       "div" : int(fields[7]),
                                       "shift" : list(map(float, fields[8:11])),
                                       "files" : fields[11:] }
                        self.phonon["animation"].append(parameters)
                    if (fields[1] == "modulations") :
                        parameters = { "q" : list(map(float, fields[2:5])),
                                       "supercell" : list(map(int,fields[5:8])), 
                                       "delta" : float(fields[8]) }
                        self.phonon["modulations"].update(parameters)

        except Exception:
            #print (line,)
            #print ("|-> line triggered exception: ") + str(Exception)
            raise
        # supercell is mandatory for all what follows
        if not self.phonon["supercell"]:
            raise Exception("no supercell specified in %s" % self.file)
        f.close()
        
