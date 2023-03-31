import numpy as np
from pykdtree.kdtree import KDTree as pyKDTree
from timeit import default_timer as timer
import sys

def outlier_filter_k(xyz, kN=50, p=1):
    '''
    Density based outlier filter. Points whose neighborhood's density is below a certain threshold are removed from the initial point cloud. Author of this code: Bodo Bookhagen
    
    Input: 
    xyz: 2D or 3D point cloud as numpy array
    kN: scalar integer that defines the neighborhood size (number of points). default: 50
    p: scalar integer that defines the percentile of the density threshold query. default: 1
    
    Output:
    xyz[idx2keep_pyKDTree]: numpy array with the filtered point cloud (points passed the density threshold)
    idx2keep_pyKDTree: numpy array, indices of points in the original point cloud that passed the density threshold query
    '''
    # array to KDTree class
    xyz_pyKDTree = pyKDTree(xyz)
    # retrieve distances to NN points
    dist_pyKDTree, _ = xyz_pyKDTree.query(xyz, k=kN)
    # dist_pyKDTree[:,-1] is point farthest away --> the smaller the distance to that point, the bigger dens_pyKDTree will be
    dens_pyKDTree =  kN / np.pi /np.power(dist_pyKDTree[:,-1],2)
    # indices for where dens_pyKDTree is over the threshold
    idx2keep_pyKDTree, = np.where(dens_pyKDTree > np.percentile(dens_pyKDTree,p))
    return xyz[idx2keep_pyKDTree], idx2keep_pyKDTree


def weight(xyz, kN=50):
    '''
    Adaptation of the density based outlier filter. Instead of removing points whose neighborhood's density is below a certain threshold, the density of each point's neighborhood is returned.
    
    Input: 
    xyz: 2D or 3D point cloud as numpy array
    kN: scalar integer that defines the neighborhood size (number of points). default: 50
    
    Output:
    dens_pyKDTree: numpy array with the density of each point's neighborhood. type: float, between 0 and 1
    '''
    # array to KDTree class
    xyz_pyKDTree = pyKDTree(xyz)
    # retrieve distances to NN points
    dist_pyKDTree, _ = xyz_pyKDTree.query(xyz, k=kN)
    # dist_pyKDTree[:,-1] is point farthest away --> the smaller the distance to that point, the bigger dens_pyKDTree will be
    dens_pyKDTree =  kN / np.pi /np.power(dist_pyKDTree[:,-1],2)
    return dens_pyKDTree

#############################################
# polyfit function

# function that handles the actual computation of polynomial fitting
def fit_calc(xy, ii, degree):
    kN_coeffs = np.empty((len(ii),degree+1))
    weight_xy = weight(xy)
    
    # linear
    if degree == 1:
        for n in np.arange(0,len(ii)):
            new_series = np.polynomial.polynomial.Polynomial.fit(xy[ii[n],0], xy[ii[n],1], deg=degree, w=weight_xy[ii[n]])
            kN_coeffs[n,:] = np.array([new_series.convert().coef[1], new_series.convert().coef[0]])
        kN_pfit_elev = kN_coeffs[:,0] * xy[:,0] + kN_coeffs[:,1]

    # quadratic    
    elif degree == 2:
        for n in np.arange(0,len(ii)):
            new_series = np.polynomial.polynomial.Polynomial.fit(xy[ii[n],0], xy[ii[n],1], deg=degree, w=weight_xy[ii[n]])
            kN_coeffs[n,:] = np.array([new_series.convert().coef[2], new_series.convert().coef[1], new_series.convert().coef[0]])
        kN_pfit_elev = kN_coeffs[:,0] * xy[:,0]**2 + kN_coeffs[:,1] * xy[:,0] + kN_coeffs[:,2]
    
    # cubic
    elif degree == 3:
        for n in np.arange(0,len(ii)):
            new_series = np.polynomial.polynomial.Polynomial.fit(xy[ii[n],0], xy[ii[n],1], deg=degree, w=weight_xy[ii[n]])
            kN_coeffs[n,:] = np.array([new_series.convert().coef[3], new_series.convert().coef[2], new_series.convert().coef[1], new_series.convert().coef[0]])
        kN_pfit_elev = kN_coeffs[:,0] * xy[:,0]**3 + kN_coeffs[:,1] * xy[:,0]**2 + kN_coeffs[:,2] * xy[:,0] + kN_coeffs[:,3]
    
    else:
        sys.exit('"degree" input is invalid, has to be an integer with the value 1, 2, or 3')
    
    return kN_coeffs, kN_pfit_elev

# function that handles the kN is array situation
def pfit_kN_wrapper(xy, kN, degree):
        xyz_pyKDTree = pyKDTree(xy) 
        # create lists that will contain coefficients and elevation
        kN_coeffs = []
        kN_pfit_elev = []
        # loop through kN array
        for g in kN:
            start_time = timer()
            global ii
            _, ii = xyz_pyKDTree.query(xy, k=g)
            kN_coeffs1, kN_pfit_elev1 = fit_calc(xy, ii, degree)
            ii = None
            kN_coeffs.append(kN_coeffs1)
            kN_pfit_elev.append(kN_pfit_elev1)
            end_timer = timer()
            print(f' kN = {g}: Time elapsed: {end_timer-start_time:.1f} seconds')
        return kN_coeffs, kN_pfit_elev 

# function sorts kN type and calls the computation function
def polyfit_numpy(xy, kN, degree):
    '''
    Function to calculate a polynomial fit with weights over a 2D point cloud. 

    Input:
    xy: 2D point cloud as numpy array. dimension: (n,2)
    kN: number of neighbors to calculate best fit for for each point. type: integer or numpy array of integers
    degree: degree of the polynomial fit. type: integer, between 1 and 3 (including 3)

    Output:
    kN_coeffs: array or list of arrays containing the coefficients of the polynomial fit for each point. the coefficients are ordered descending by degree. dimension of array(s): (length of xy, degree+1) 
    kN_pfit_elev: array or list of arrays containing the function value (--> elevation) of the polynomial fit for each point. dimension of array(s): (length of xy, 1)
    '''
    # check if kN is single integer or array
    if type(kN) == int:
        xyz_pyKDTree = pyKDTree(xy)
        _, ii = xyz_pyKDTree.query(xy, k=kN)
        # calculate coefficients and fitted elevation
        kN_coeffs, kN_pfit_elev = fit_calc(xy, ii, degree)
    elif type(kN) == np.ndarray:
        # catch if xy is shorter than biggest neighborhood, prevents kernel from dying
        if len(xy) < kN[-1]:
            sys.exit('xy too short')
        # calculate coefficients and fitted elevation
        kN_coeffs, kN_pfit_elev = pfit_kN_wrapper(xy, kN, degree)
    # if not, stop function
    else:
        sys.exit('unsupported kN type (has to be integer or np.array)')
    # final return: coefficients and fitted elevation as either np.array or list of np.arrays
    return kN_coeffs, kN_pfit_elev