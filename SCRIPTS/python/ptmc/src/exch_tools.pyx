import sys
import numpy as np
cimport numpy as np
cimport cython

#cimport _random_displace
cdef extern:
    int calculate_diffusion(long int * paths, long int *diffusion, long int *rdiff,
        long int nreplicas, long int ntimes, long int *maxelapsed, long int *maxelapsed_cum)

def calculate_diffusion_cython(np.ndarray[long int, ndim=2, mode="c"] paths
                ):
    ntimes = len(paths[:,0])
    nreplicas = len(paths[0,:])
    cdef np.ndarray[long int, ndim=1, mode="c"] diffusion = np.zeros([nreplicas*nreplicas], dtype=np.integer)
    cdef np.ndarray[long int, ndim=1, mode="c"] rdiff = np.zeros(nreplicas*nreplicas, dtype=np.integer)
    cdef np.ndarray[long int, ndim=1, mode="c"] maxelapsed = np.zeros(1, dtype=np.integer)
    cdef np.ndarray[long int, ndim=1, mode="c"] maxelapsed_cum = np.zeros(1, dtype=np.integer)
    sys.stdout.flush()
    cdef np.ndarray[long int, ndim=1, mode="c"] paths1d = paths.flatten()

    sys.stdout.flush()
    calculate_diffusion(
            <long int*>paths1d.data, 
            <long int*>diffusion.data, 
            <long int*>rdiff.data, 
            nreplicas, ntimes,
            <long int*>maxelapsed.data, 
            <long int*>maxelapsed_cum.data
                      )
    return diffusion.reshape([nreplicas,nreplicas]), rdiff.reshape([nreplicas, nreplicas]), maxelapsed[0], maxelapsed_cum[0]


@cython.boundscheck(False) # turn of bounds-checking for entire function
def calculate_occupation_probability_cython(np.ndarray[long int, ndim=2, mode="c"] paths):
    cdef unsigned int ntimes = len(paths[:,0])
    cdef unsigned int nreplicas = len(paths[0,:])
    cdef np.ndarray[double, ndim=2, mode="c"] occupation_prob = np.zeros([nreplicas, nreplicas])
    cdef unsigned int i, r, p, t
    
    for t in range(ntimes):
        for r in xrange(nreplicas):
            i = paths[t,r]
            occupation_prob[r, i] += 1.
    
    return occupation_prob

@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.wraparound(False) # turn of bounds-checking for entire function
def path_from_permutations(np.ndarray[long int, ndim=2, mode="c"] permutations):
    cdef unsigned int nperm = len(permutations[:,0])
    cdef unsigned int nreplicas = len(permutations[0,:])
    cdef np.ndarray[long int, ndim=2, mode="c"] paths = np.zeros(np.shape(permutations), dtype=np.integer)
    cdef np.ndarray[long int, ndim=2, mode="c"] replicas = np.zeros(np.shape(permutations), dtype=np.integer)
    cdef np.ndarray[long int, ndim=1, mode="c"] exchanges_counts = np.zeros(nreplicas, dtype=np.integer)
    
    paths[0,:] = range(nreplicas)
    replicas[0,:] = range(nreplicas)

    
    cdef unsigned int t, r, old_pos, pos, rep
#    cdef int t

    for t in range(1,nperm):
        for r in range(nreplicas): 
            old_pos = permutations[t-1,r]
            replicas[t,r] = replicas[t-1, old_pos]
        for pos in range(nreplicas):
            rep = replicas[t,pos]
            paths[t,rep] = pos
        for r in range(nreplicas):
            if permutations[t,r] != r:
                exchanges_counts[r] += 1
    
    return paths, replicas, exchanges_counts
    
