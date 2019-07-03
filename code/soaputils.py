import numpy as np
import soaplite
import genBasis
import ase
import matplotlib.pyplot as p
from numpy.linalg import norm, svd
from ase.visualize import view
import scipy.spatial.distance as sd


def rand_pos(atoms_in, seed=None): # function to randomize atom positions. New positions will be inside Unit cell
    atoms = atoms_in.copy()
    pos = atoms.get_positions()
    shape = pos.shape
    np.random.seed(seed)
    ran_pos = np.random.random_sample(shape)
    cell = atoms.get_cell()
    for i in np.arange(shape[0]):
        pos[i] = np.matmul(cell,ran_pos[i])
    atoms.set_positions(pos)
    return atoms

def limit_pos(atoms_obj): # funtion to delete atoms outside of unit cell, only works with orthorhombic cells (for now)
    atoms = atoms_obj.copy()
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

def soap_norm(pos, atoms_obj, myAlphas=0, myBetas=0, rCut=10.0, NradBas=5, Lmax=5, pbc=False, maximize=False, showProgress=False):
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
    if maximize:
        return -norm(mat)
    else:
        return norm(mat)
    
def soap_norm_LJ(pos, atoms_obj, sigma=2, eps=1 ,myAlphas=0, myBetas=0, rCut=10.0, NradBas=5, Lmax=5, pbc=False, maximize=False, showProgress=False):
    # calculates the matrix-norm of the SOAP-matrix from a (flattened) positions-array and an atoms object.
    # SOAP basisfunctions have to be calculated beforehand and parsed to the function, rCut, NradBas and Lmax may also be parsed
    if myAlphas.all() == 0 or myBetas.all() == 0:
        myAlphas, myBetas = genBasis.getBasisFunc(rCut, NradBas)
    pos_ini = atoms_obj.get_positions()
    pos = np.reshape(pos, pos_ini.shape)
    atoms_obj.set_positions(pos)
    cost = cost_LJ_trunc(atoms_obj, sigma, eps, rCut, pbc)
    if pbc:
        mat = soaplite.get_periodic_soap_structure(atoms_obj, myAlphas, myBetas, rCut, NradBas, Lmax)
    else:
        mat = soaplite.get_soap_structure(atoms_obj, myAlphas, myBetas, rCut, NradBas, Lmax)
    if showProgress:
        print(norm(mat)+cost) # can be used to show progress, but slows down function calls somewhat
    if maximize:
        return -norm(mat)+cost
    else:
        return norm(mat)+cost
        
def svd_lp(pos, atoms_obj, order=2, myAlphas=0, myBetas=0, rCut=10.0, NradBas=5, Lmax=5, pbc=False, showProgress=False):
    # calculates and returns norm of singular-value-vector of SOAP-matrix, works similar to soap_norm
    # the order of the norm is 2 by default, but can be set in variable 'order'
    #if myAlphas.all() == 0 or myBetas.all() == 0:
    if myAlphas.all() or myBetas.all():
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
        print(norm(s, ord=order)) # can be used to show progress, but slows down function calls somewhat
    return  norm(s, ord=order)

def norm_block(pos, atoms_obj, order=2, frac=0.3, myAlphas=0, myBetas=0, rCut=10.0, NradBas=5, Lmax=5, pbc=False, showProgress=False):
    # SOAP-matrix is divided into blocks depending on combination of atomic number. The norm of the last rows of the blocks is calculated (fraction of rows is put in with frac, so last 2/3 of rows would be frac = 1/3), rest works similar to soap_norm
    # the order of the norm is 2 by default, but can be set in variable 'order'
    #if myAlphas.all() == 0 or myBetas.all() == 0:
    if myAlphas.all() or myBetas.all():
        myAlphas, myBetas = genBasis.getBasisFunc(rCut, NradBas)
    pos_ini = atoms_obj.get_positions()
    pos = np.reshape(pos, pos_ini.shape)
    atoms_obj.set_positions(pos)
    if pbc:
        mat = soaplite.get_periodic_soap_structure(atoms_obj, myAlphas, myBetas, rCut, NradBas, Lmax)
    else:
        mat = soaplite.get_soap_structure(atoms_obj, myAlphas, myBetas, rCut, NradBas, Lmax)
        
    symbols = atoms_obj.get_chemical_symbols() # here the blocks are seperated and the 'block norm' is calculated
    uSym = len(np.unique(symbols))
    numBlocks = int(uSym*(uSym + 1)/2)
    rowsBlock = int(np.shape(mat)[1]/numBlocks)
    mat0 = np.zeros((np.shape(mat)[0],1))
    for i in np.arange(numBlocks):
        mat0 = np.append(mat0,mat[:,i*rowsBlock+int(frac*rowsBlock):(i+1)*rowsBlock],axis=1)
    
    if showProgress:    
        print(norm(mat0, ord=order)) # can be used to show progress, but slows down function calls somewhat
    return  norm(mat0, ord=order)

def show_res(atoms_obj, pos, sigma=1, eps=1, myAlphas=0, myBetas=0, rCut=10.0, NradBas=5, Lmax=5, pbc=False):
    # shows the result of the minimization in form of the view from ase
    #if myAlphas.all() == 0 or myBetas.all() == 0:
    if myAlphas.all() or myBetas.all():
        myAlphas, myBetas = genBasis.getBasisFunc(rCut, NradBas)
    atoms = atoms_obj.copy()
    pos_ini = atoms.get_positions()
    pos = np.reshape(pos, pos_ini.shape)
    atoms.set_positions(pos)
    if pbc:
        mat = soaplite.get_periodic_soap_structure(atoms, myAlphas, myBetas, rCut, NradBas, Lmax)
    else:
        mat = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut, NradBas, Lmax)
    s = svd(mat.transpose(), full_matrices=False, compute_uv=False)
    norm_b = norm_block(pos, atoms, 2, 0.3, myAlphas, myBetas, rCut, NradBas, Lmax, pbc)
    #print('Matrix norm: %f' %norm(mat))
    N = len(atoms.get_positions())*3
    print('Number of Atoms: %i' %(N/3))
    dist = atoms.get_all_distances(mic=pbc)
    np.fill_diagonal(dist,10.0)
    amin = np.amin(dist)
    cost = cost_LJ_trunc(atoms, sigma, eps, rCut, pbc)
    print('Minimal pair-distance: %f' %amin)
    print('Singular-Value norm: %f' %norm(s, ord=1))
    print('Full Matrix norm: %f' %norm(mat, ord=2))
    print('Block matrix norm: %f' %norm_b)
    print('LJ-Cost: %f' %cost)
    #p.matshow(mat)
    #p.semilogy(s)
    f, (ax1, ax2) = p.subplots(2, 1)
    ax1.matshow(mat)
    ax2.semilogy(s)
    view(atoms)
    return(norm(mat, ord=2), cost)
    

    
def lim_overlap_ran(atoms_obj, dmin, max_it=10000,seed=None):
    # tries to make atoms in the cell not overlap, by assigning new random position of atom inside cell. Minimum distance between atompairs is dmin
    obj = atoms_obj.copy()
    overlap = True
    cell = obj.get_cell()
    it = 0
    np.random.seed(seed)
    #N = len(obj.get_atomic_numbers())
    while(overlap and it < max_it):
        dist = obj.get_all_distances(mic=True)
        #print(obj)
        np.fill_diagonal(dist,1000.0)
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
    if it >= 1000:
        print('Maximum number of iterations (%i) to limit overlap exceeded' %max_it)
    return obj

def lim_overlap(atoms_obj, dmin, max_it=10000,seed=None):
    # makes atoms not overlap by setting the distance between the closest atoms and iterating on that, works better than the random algorithm but only up to 3.3A
    obj = atoms_obj.copy()
    overlap = True
    cell = obj.get_cell()
    it = 0
    np.random.seed(seed)
    #N = len(obj.get_atomic_numbers())
    while(overlap and it < max_it):
        dist = obj.get_all_distances(mic=True)
        #print(obj)
        np.fill_diagonal(dist,1000.0)
        #print(np.amin(dist))
        if np.amin(dist) > dmin:
            overlap = False
        else:
            index = np.unravel_index(np.argmin(dist), np.shape(dist))
            pos = obj.get_positions()
            shape = pos.shape
            obj.set_distance(index[0], index[1], dmin, mic=True)
            it = it + 1
    if it >= 1000:
        print('Maximum number of iterations (%i) to limit overlap exceeded' %max_it)
    return obj

def cost_LJ(atoms, sigma, eps):
    arr = atoms.get_positions()
    r = sd.pdist(arr)/sigma
    np.power(r, -6, out=r)
    return np.sum(r**2 - r)*4*eps
    
def cost_LJ_trunc(atoms, sigma=1, eps=1, rCut=5, pbc=False):
    if pbc:
        r = atoms.get_all_distances(mic=True)
        r.flatten()
        r = r[(r != 0.0)]
    else:
        arr = atoms.get_positions()
        r = sd.pdist(arr)
    shift = np.power(rCut, -6)
    r = r[(r < rCut)]
    np.power(r/sigma, -6, out=r)
    return np.sum(r**2 - r - shift**2 + shift)*4*eps
    
    
    
    
    
    
    
    
    
    
    