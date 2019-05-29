import numpy as np
import soaplite
import genBasis
import ase
import matplotlib.pyplot as p
from numpy.linalg import norm, svd
from ase.visualize import view


def rand_pos(atoms_in): # function to randomize atom positions. New positions will be inside Unit cell
    atoms = atoms_in.copy()
    pos = atoms.get_positions()
    shape = pos.shape
    ran_pos = np.random.random_sample(shape)
    cell = atoms.get_cell()
    for i in np.arange(shape[0]):
        pos[i] = np.matmul(cell,ran_pos[i])
    atoms.set_positions(pos)
    return atoms

def limit_pos(atoms): # funtion to delete atoms outside of unit cell, only works with orthorhombic cells (for now)
    cell = atoms.get_cell()
    pos = atoms.get_positions()
    N = np.shape(pos)[0]
    x_max = cell[0,0]
    y_max = cell[1,1]
    z_max = cell[2,2]
    i = 0
    while(i < N-1):
        pos = atoms.get_positions()
        N = np.shape(pos)[0]
        if pos[i,0] > cell[0,0] or pos[i,1] > cell[1,1] or pos[i,2] > cell[2,2]:
            atoms.pop(i)
        else:
            i = i + 1
    return atoms

def soap_norm(pos, atoms_obj, myAlphas=0, myBetas=0, rCut=10.0, NradBas=5, Lmax=5, pbc=False, showProgress=False):
    # calculates the matrix-norm of the SOAP-matrix from a (flattened) positions-array and an atoms object.
    # SOAP basisfunctions have to be calculated beforehand and parsed to the function, rCut, NradBas and Lmax may also be parsed
    if myAlphas.all() == 0 or myBetas.all() == 0:
        myAlphas, myBetas = genBasis.getBasisFunc(rCut, NradBas)
    pos_ini = atoms_obj.get_positions()
    pos = np.reshape(pos, pos_ini.shape)
    atoms_obj.set_positions(pos)
    if pbc:
        mat = soaplite.get_periodic_soap_structure(atoms_obj, myAlphas, myBetas, rCut, NradBas, Lmax)
    else:
        mat = soaplite.get_soap_structure(atoms_obj, myAlphas, myBetas, rCut, NradBas, Lmax)
    if showProgress:    
        print(norm(mat)) # can be used to show progress, but slows down function calls somewhat
    return norm(mat)

def svd_l2(pos, atoms_obj, myAlphas=0, myBetas=0, rCut=10.0, NradBas=5, Lmax=5, pbc=False, showProgress=False):
    # calculates and returns norm of singular-value-vector of SOAP-matrix, works similar to soap_norm
    if myAlphas.all() == 0 or myBetas.all() == 0:
        myAlphas, myBetas = genBasis.getBasisFunc(rCut, NradBas)
    pos_ini = atoms_obj.get_positions()
    pos = np.reshape(pos, pos_ini.shape)
    atoms_obj.set_positions(pos)
    if pbc:
        mat = soaplite.get_periodic_soap_structure(atoms_obj, myAlphas, myBetas, rCut, NradBas, Lmax)
    else:
        mat = soaplite.get_soap_structure(atoms_obj, myAlphas, myBetas, rCut, NradBas, Lmax)
    s = svd(mat.transpose(), full_matrices=False, compute_uv=False)
    if showProgress:    
        print(norm(s, ord=1)) # can be used to show progress, but slows down function calls somewhat
    return  norm(s, ord=1)

def show_res(atoms_obj, pos, myAlphas=0, myBetas=0, rCut=10.0, NradBas=5, Lmax=5, pbc=False):
    # shows the result of the minimization in form of the view from ase
    if myAlphas.all() == 0 or myBetas.all() == 0:
        myAlphas, myBetas = genBasis.getBasisFunc(rCut, NradBas)
    atoms = atoms_obj.copy()
    pos_ini = atoms.get_positions()
    pos = np.reshape(pos, pos_ini.shape)
    atoms.set_positions(pos)
    if pbc:
        mat = soaplite.get_periodic_soap_structure(atoms_obj, myAlphas, myBetas, rCut, NradBas, Lmax)
    else:
        mat = soaplite.get_soap_structure(atoms_obj, myAlphas, myBetas, rCut, NradBas, Lmax)
    s = svd(mat.transpose(), full_matrices=False, compute_uv=False)
    #print('Matrix norm: %f' %norm(mat))
    print('Singular-Value norm: %f' %norm(s, ord=1))
    #p.matshow(mat)
    #p.semilogy(s)
    f, (ax1, ax2) = p.subplots(2, 1)
    ax1.matshow(mat)
    ax2.semilogy(s)
    view(atoms)
    
def lim_overlap(atoms_obj, dmin=1.1):
    # tries to make atoms in the cell not overlap, by assigning new random position of atom inside cell. Minimum distance between atompairs is dmin
    obj = atoms_obj.copy()
    overlap = True
    cell = obj.get_cell()
    it = 0
    #N = len(obj.get_atomic_numbers())
    while(overlap and it < 1000):
        dist = obj.get_all_distances(mic=True)
        #print(obj)
        np.fill_diagonal(dist,10.0)
        #print(np.amin(dist))
        if np.amin(dist) > dmin:
            overlap = False
        else:
            index = np.unravel_index(np.argmin(dist), np.shape(dist))[0]
            pos = obj.get_positions()
            shape = pos.shape
            ran_pos = np.random.random_sample(shape)
            pos[index] = np.matmul(cell,ran_pos[index])
            obj.set_positions(pos)
            it = it + 1
    return obj

    