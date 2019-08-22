import ase
import numpy as np
from ase.neighborlist import neighbor_list



def computePPDF(struct,binSize=0.3,numBins=20, projectedAxis=None):
    '''The function computePPDF will compute the partial pair distribution function for the complete structure. 
    '''
    # First, we want to make sure that the cell has periodic boundary conditions.
    unit_cell = struct.get_cell() 
    pbc = struct.get_pbc()
    useMIC = 1
    if ~np.all(pbc):
        print("This structure does not have all boundaries being periodic (%d %d %d)"%(pbc[0],pbc[1],pbc[2]))
    
    # Next, we must compute the distance of every atom to each other atom, because we then want to generate a 
    # histogram of those. The function get_all_distances() returns a matrix of the distances between each atom 
    # every other atom. The parameter 'mic' specifies wether to use the minimum image convention (MIC). 
    # If we set it to true, then priodic boundary conditions are effectively being used, this means the distance 
    # of an atom to its nearest image in any of the surrounding identical copies of the simulation box
    # is being computed.
    if (projectedAxis == None):
        dist = struct.get_all_distances(mic=useMIC)
    else:
        R0 = struct.get_array('positions')
        R = np.copy(R0)
        R[:,projectedAxis] = 0
        struct.set_array('positions', R)
        dist = struct.get_all_distances(mic=useMIC)
        # copy the original atom positions back into their place:
        struct.set_array('positions', R0)
         
    # In order to not include the distance between each atom with itself, we will set the diagonal elements
    # of this matrix to -1:
    np.fill_diagonal(dist, -1.0)
    
    # Compute a histogram of values in dist, but sorted according to the Z-numbers (in reverse order).
    sym = struct.get_chemical_symbols()
    Natom = len(sym)  # number of atoms = number of symbols
    uSym,uZ = extractUniqueSymbols(struct)
    index = np.zeros((Natom,))
    for j in range(Natom):
        index[j] = uSym.index(sym[j])
    # The array 'index' now specifies the position of each chemical symbol in the list of unique symbols.
    # For example, if the list of symbols is as follows:
    # ['Cu', 'Cu', 'Zr', 'Zr', 'Cu', 'Zr', 'Zr', 'Zr', 'Zr', 'Cu', 'Zr', 'Cu', 'Cu', ...
    # the index array will hold the following entries:
    # [ 1.  1.  0.  0.  1.  0.  0.  0.  0.  1.  0.  1.  1.  ...
    # The highest Z-number has the lowest index in this array, because it has the strongest contribution
    # to the potential
    
    # Next, we will prepare a matrix of the column index in which 
    NuSym = len(uSym)  # number of unique elements
    indexArray = np.zeros((Natom,Natom,2))
    indexArray[:,:,0] = np.tile(index.reshape(1, -1),[Natom,1])
    indexArray[:,:,1] = np.tile(index.reshape(-1, 1),[1,Natom])
    indexArray = np.sort(indexArray,axis=2)
    # When creating the index array we will put the Zr-Cu distances in the same histogram as the 
    # Cu-Zr distances. For this reason we need to subtract from the index of the array those possibilities
    # that account for the opposite sequence:
    symArray = np.zeros((Natom,Natom),dtype=int)
    symArray[:,:] = NuSym*indexArray[:,:,0]+indexArray[:,:,1]-(indexArray[:,:,0]+1)*indexArray[:,:,0]/2
    #print(symArray[0:5,0:5])
    
    # Compute the number of different histograms we need:
    numHist = np.max(symArray,axis=None)+1
    
    histRange = (0,numBins*binSize)
    hist = np.zeros((numBins,numHist))
    
    for combination in range(numHist):
        h, hEdge = np.histogram(np.extract(symArray == combination, dist),bins=numBins,range=histRange)
        hist[:,combination] = h
        
    # Now, we also want to normalize the positions to the size of the unit cell along the dimension we have 
    # integrated over:
    if (projectedAxis != None):
        hist /= unit_cell[projectedAxis,projectedAxis]
    distBins = 0.5*(hEdge[1:]+hEdge[0:-1])

    return hist, distBins

def computeRDF(struct,binSize=0.3,numBins=20):
    hist, distBins = computePPDF(struct,binSize,numBins)
    rho = density(struct)
    dr = distBins[1]-distBins[0]
    rhist = [hist[i]/(rho*4*np.pi*distBins[i]**2*dr) for i in np.arange(len(distBins))]
    return rhist, distBins

def simpleRDF(struct,numBins=100,rCut=10):
    d = neighbor_list('d', struct, rCut)
    h, bin_edges = np.histogram(d, bins=numBins)
    rdf = h/(4*np.pi/3*(bin_edges[1:]**3 - bin_edges[:-1]**3)) * struct.get_volume()/len(struct)**2
    distBins = 0.5*(bin_edges[1:]+bin_edges[:-1])
    return rdf, distBins
    

def extractUniqueSymbols(struct):    
    sym = struct.get_chemical_symbols()
    Z = struct.get_atomic_numbers()
    uniqueSymbols = list(set(sym))
    uniqueZ = list(set(Z))
    uniqueEl = sorted(uniqueSymbols, key=lambda uniqueZ: uniqueZ, reverse=True)
    # Now we want to sort the symbols coording to the Z-numbers in descending order
    uniqueZ.sort()
    
    return uniqueEl, uniqueZ

def density(struct):
    cell = struct.get_cell()
    x = cell[:,0]
    y = cell[:,1]
    z = cell[:,2]
    vol = np.dot(np.cross(x,y),z)
    n = len(struct.get_chemical_symbols())
    return n/vol
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    