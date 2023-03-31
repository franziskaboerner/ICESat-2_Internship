import numpy as np
from pykdtree.kdtree import KDTree as pyKDTree
from numba import njit, prange
import sys

# numba eigenvalue function without kN loop

def numba_wrapper_eva(xy, ii):
    @njit(parallel=True)
    def numba_eva(xy, ii):
        # prange = parallel range, run loops in parallel, for as long as input is long
        # empty array that will store each point's mean eigenvalues
        # dimensions: row amount of input ROWS and 2 COLUMNS
        ev_all = np.empty((len(ii),2))
        std_all = np.empty((len(ii),2))
        for n in prange(len(ii)):
            # create array that contains the mean of each dimension of the NN indexed points 
            # dimension of c: (1,2)
            c = np.array([xy[ii[n], 0].mean(), xy[ii[n], 1].mean()])
            # subtract each value by the corresponding mean --> shifted so center of mass is at (0,0)
            q = xy[ii[n], :] - c
            # define number of iterations for randomizer
            rand_it = 50
            # empty array to store eigenvalues from the randomizer loop
            # is newly set up for every point n
            ev = np.empty((rand_it,2))
            for j in prange(rand_it):
                # randomly choose ~75% of points to calculate covariance matrix from 
                # --> get indices to subset both distance and height
                indices = np.random.choice(q[:,1].shape[0], int(np.floor(0.75*q.shape[0])), replace=False)
                # index the distance and height
                hindex = q[:,1][indices]
                distindex = q[:,0][indices]
                # put sampled dist and height back in one matrix (already in 2,n shape --> does not need to be transposed)
                q_rand = np.row_stack((distindex, hindex))
                # calculate eigenvalues from covariance matrix
                # dim: (1,2)
                evalue = np.linalg.eigh(np.cov(q_rand))[0]
                # write eigenvalues to ev array
                # this array includes all j-range EV pairs for the j'th point in the prange loop
                # dim: (j-range, 2)
                ev[j,:] = evalue
            # calculate the mean EVs and StDev for point n and its neighbors, add to arrays for each point
            ev_all[n,:] = np.array([np.mean(ev[:,0]), np.mean(ev[:,1])])
            std_all[n,:] = np.array([np.std(ev[:,0]), np.std(ev[:,1])])
        # return the main EV and StDev array
        return ev_all, std_all
    return numba_eva(xy, ii)

def kN_wrapper(xy):
        # put xyz array in KDTree class 
        xyz_pyKDTree = pyKDTree(xy)
        # empty lists for eigenvalues and standard deviation    
        kN_list_ev = []
        kN_list_std = []
        for g in kN:
            global ii
            # query returns "closest_dists_res, closest_idxs_res"
            # closest_dists_res = cartesian/euclidian distances
            # closest_idxs_res = indices of NN points
            # only take point indices of each point's kN nearest 
            # neighbors
            _, ii = xyz_pyKDTree.query(xy, k=g)
            # define n as number of rows of xyz array
            n = xy.shape[0]
            # call the numba eigenvalue wrapper to calculate eigenvalues
            # of all chosen NN points
            ev_all, std_all = numba_wrapper_eva(xy, ii)
            # clear ii
            ii = None
            # append ev/std of current kN to lists
            kN_list_ev.append(ev_all)
            kN_list_std.append(std_all)
        return kN_list_ev, kN_list_std 

def evalues_pyKDTree_kN_loop(xy):
    '''
    Function to calculate eigenvalues of 2D point cloud arrays based on point nearest neighborhood (10, 50, 100, 150, 250, 400 neighbors). For each point, eigenvalue calculation is performed 50 times with a random sample of 75% of the point's neighbors.   

    Input:
    xy: 2D point cloud as numpy array. dimension: (n,2)

    Output: 
    kN_list_ev: list (length: 6) of numpy arrays with eigenvalues of each point, based on eigenvalues of the point's neighborhood
    kN_list_ev: list (length: 6) of numpy arrays with standard deviation of eigenvalues of each point from the randomizer
    '''

    # initiate neighborhood sizes
    global kN   
    kN = np.array([10, 50, 100, 150, 250, 400])
    # catch if xy is shorter than biggest neighborhood, prevents kernel from dying
    if len(xy) < kN[-1]:
        sys.exit('xy too short')
    # for loop via numba to retrieve eigenvalues and standard deviation
    kN_list_ev, kN_list_std = kN_wrapper(xy)
    # return the eigenvalues and stdev
    return kN_list_ev, kN_list_std

