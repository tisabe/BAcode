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
This file computes a batch of optimizations with different numbers of atoms in parallel. This is done to test the scalability with size of the structure.
'''


rCut = 10
NradBas = 5
Lmax = 5
myAlphas, myBetas = genBasis.getBasisFunc(rCut, NradBas)
    
N_list = [10,20,30,40,50,60,70,80,90,100,150,200]


def calc(i):
    Natoms = i
    atoms = su.gen_struct(Natoms, seed=50)
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
    filename = "batchC"+str(i)+".cfg"
    f = open("res_structs/size_test/log.txt","a")
    f.write("N: " + str(i) + " nit: " + str(nit) + " nfev: " + str(nfev) + " t: " + str(dt) + " msg: " + str(msg) + "\n")
    f.close()
    ase.io.write("res_structs/size_test/" + filename, atoms_res)
    return 0 # return value can be ignored

pool = mp.Pool(processes=6)
results = [pool.apply_async(calc, args=(i,)) for i in N_list]; output = [p.get() for p in results]
pool.terminate()

