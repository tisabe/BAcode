import numpy as np
import soaplite
import genBasis
import ase
import matplotlib.pyplot as p
from numpy.linalg import norm, svd
import soaputils as su
import scipy.optimize as op
import multiprocessing as mp
import time

'''
This file computes a batch of optimizations with different SOAP-parameters in parallel. This is done to test effect of different parameters on optimisation result.
'''
l_list = [0,1,3,5,7,9,5,5,5,5,5,5,5,5,5,5,5,5]
n_list = [5,5,5,5,5,5,2,4,5,8,10,12,5,5,5,5,5,5]
r_list = [10,10,10,10,10,10,10,10,10,10,10,10,2,4,5,8,10,12]

Natoms = 70
atoms = su.gen_struct(Natoms, seed=50)

def calc70(i):
    rCut = r_list[i]
    NradBas = n_list[i]
    Lmax = l_list[i]
    myAlphas, myBetas = genBasis.getBasisFunc(rCut, NradBas)
    
    cell_len = atoms.get_cell()[0,0]
    x0 = atoms.get_positions()
    bounds_obj = [(0,cell_len)]*Natoms*3
    op_options = {'maxiter': 50000, 'disp': True}
    t0 = time.time()
    
    res_obj = op.minimize(su.svd_norm, x0, method='L-BFGS-B',args=(atoms, myAlphas, myBetas, rCut, NradBas, Lmax, True), bounds=bounds_obj,  options=op_options)
    
    t1 = time.time()
    dt = t1 - t0
    xopt = res_obj.x
    nit = res_obj.nit
    nfev = res_obj.nfev
    msg = res_obj.message
    atoms_res = atoms.copy()
    atoms_res.set_positions(np.reshape(xopt,(-1,3)))
    
    # the resulting struct is saved in a folder
    filename = "batch" + "l" + str(l_list[i]) + "n" + str(n_list[i]) + "r" + str(r_list[i]) + ".cfg"
    f = open("res_structs/param_test/msg1.txt","a")
    f.write("run " + str(i) + " nit: " + str(nit) + " nfev: " + str(nfev) + " t: " + str(dt) + " msg: " + str(msg) + "\n")
    f.close()
    ase.io.write("res_structs/param_test/" + filename, atoms_res)
    return 0 # return value can be ignored

indices = np.arange(len(l_list))
pool = mp.Pool(processes=6)
results = [pool.apply_async(calc70, args=(i,)) for i in indices]; output = [p.get() for p in results]
pool.terminate()










