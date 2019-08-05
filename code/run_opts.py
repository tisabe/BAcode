import numpy as np
import soaplite
import genBasis
import ase
import matplotlib.pyplot as p
from numpy.linalg import norm, svd
import soaputils as su
import scipy.optimize as op
import multiprocessing as mp

'''
This file computes a batch of optimizations with different starting positions in parallel. This is done to show reproducibility between optimizations with different starting points.
'''
# parameters for soap-calculation are set
rCut = 10.0
NradBas = 5
Lmax = 5
myAlphas, myBetas = genBasis.getBasisFunc(rCut, NradBas)

# this function performs an optimization of the normalised svd with 70 atoms. The parameter i is a seed-offset to produce different starting values for the random structure generation
def calc70(i):
    Natoms = 70
    atoms = su.gen_struct(Natoms, seed=50+i)
    cell_len = atoms.get_cell()[0,0]
    x0 = atoms.get_positions()
    bounds_obj = [(0,cell_len)]*Natoms*3
    op_options = {'maxiter': 50000, 'disp': True}
    res_obj = op.minimize(su.svd_norm, x0, method='L-BFGS-B',args=(atoms, myAlphas, myBetas, rCut, NradBas, Lmax, True), bounds=bounds_obj,  options=op_options)
    xopt = res_obj.x
    atoms_res = atoms.copy()
    atoms_res.set_positions(np.reshape(xopt,(-1,3)))
    # the resulting struct is saved in a folder
    filename = "batch" + str(i) + ".cfg"
    ase.io.write("res_structs/x0_test/" + filename, atoms_res)
    return 0 # return value can be ignored

# here the multiprocessing pool is initialised and the calculation started
seeds = np.arange(10)
pool = mp.Pool(processes=6)
results = [pool.apply_async(calc70, args=(i,)) for i in seeds]; output = [p.get() for p in results]
pool.terminate()






