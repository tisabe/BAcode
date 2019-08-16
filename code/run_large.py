import numpy as np
import soaplite
import genBasis
import ase
import soaputils as su
import scipy.optimize as op
import time

Natoms = 100
atoms = su.gen_struct(Natoms, seed=50)

rCut = 10
NradBas = 5
Lmax = 5
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
f = open("res_structs/log_run_large.txt","a")
f.write(" nit: " + str(nit) + " nfev: " + str(nfev) + " t: " + str(dt) + " msg: " + str(msg) + "\n")
f.close()
ase.io.write("res_structs/C100.cfg", atoms_res)
