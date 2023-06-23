import sys
import random
import numpy as np
import scanpy as sc
import paste as pst
import pandas as pd

slide1_counts = "data/tsv/raw-counts-sc1.tsv"
slide1_coords = "data/tsv/coords-slide1.tsv"
slide1_raw_coords = "data/tsv/raw-coords-slide1.tsv" 
slide2_counts = "data/tsv/raw-counts-sc2.tsv"
slide2_coords = "data/tsv/coords-slide2.tsv"
slide2_raw_coords = "data/tsv/raw-coords-slide2.tsv"

def spatial_scanpy(counts, coords):
    "Returns a scanpy spatial object. Requires raw and already filtered counts."
    scanpyobj = sc.read_csv(counts, delimiter = "\t", first_column_names = True)
    scanpyobj = scanpyobj.transpose()
    scanpyobj.obsm["spatial"] = np.genfromtxt(coords, delimiter = "\t", \
                                              skip_header = 1, usecols = (1, 2))
    return scanpyobj

def generalized_procrustes_analysis(X, Y, pi, output_params = False, matrix = False):
    """
    From https://github.com/raphael-group/paste visualization.py
    Finds and applies optimal rotation between spatial coordinates of two layers (may also do a reflection).
    Args:
        X: np array of spatial coordinates (ex: sliceA.obs['spatial'])
        Y: np array of spatial coordinates (ex: sliceB.obs['spatial'])
        pi: mapping between the two layers output by PASTE
        output_params: Boolean of whether to return rotation angle and translations along with spatial coordiantes.
        matrix: Boolean of whether to return the rotation as a matrix or an angle.
    Returns:
        Aligned spatial coordinates of X, Y, rotation angle, translation of X, translation of Y.
    """
    assert X.shape[1] == 2 and Y.shape[1] == 2

    tX = pi.sum(axis=1).dot(X)
    tY = pi.sum(axis=0).dot(Y)
    X = X - tX
    Y = Y - tY
    H = Y.T.dot(pi.T.dot(X))
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T.dot(U.T)
    Y = R.dot(Y.T).T
    if output_params and not matrix:
        M = np.array([[0,-1],[1,0]])
        theta = np.arctan(np.trace(M.dot(H))/np.trace(H))
        return X,Y,theta,tX,tY
    elif output_params and matrix:
        return X, Y, R, tX, tY
    else:
        return X,Y
    

## CODE ##
# Random seed
random.seed(1)

# Create Scanpy objects
slide1 = spatial_scanpy(slide1_counts, slide1_coords)
slide2 = spatial_scanpy(slide2_counts, slide2_coords)
slide1_raw = spatial_scanpy(slide1_counts, slide1_raw_coords)
slide2_raw = spatial_scanpy(slide2_counts, slide2_raw_coords)

ref = slide1
target = slide2
ref_raw = slide1_raw
target_raw = slide2_raw

# Pair-wise alignment of slides
pis = pst.pairwise_align(target, ref)
pis_raw = pst.pairwise_align(target_raw, ref_raw)

# Compute aligned coordinates
Target, Ref, angle, tTarget, tRef = generalized_procrustes_analysis(target.obsm["spatial"], ref.obsm["spatial"], pis, \
                                                                    output_params = True, matrix = True)
Ref = Ref.dot(np.linalg.inv(angle.T))

Target_raw, Ref_raw, angle_raw, tTarget_raw, tRef_raw = generalized_procrustes_analysis(target_raw.obsm["spatial"], ref_raw.obsm["spatial"], pis_raw, \
                                                                    output_params = True, matrix = True)
Ref_raw = Ref_raw.dot(np.linalg.inv(angle_raw.T))

# Craete pandas dataframes
colnames = ["adj_x", "adj_y"]
target_df = pd.DataFrame(Target, index = target.obs_names, columns = colnames)
ref_df = pd.DataFrame(Ref, index = ref.obs_names, columns = colnames)

target_raw_df = pd.DataFrame(Target_raw, index = target_raw.obs_names, columns = colnames)
ref_raw_df = pd.DataFrame(Ref_raw, index = ref_raw.obs_names, columns = colnames)

# Save dataframes
ref_df.to_csv("data/tsv/aligned_coords_sc1.tsv", sep = "\t", header = True, index = True)
target_df.to_csv("data/tsv/aligned_coords_sc2.tsv", sep = "\t", header = True, index = True)

ref_raw_df.to_csv("data/tsv/aligned_raw_coords_sc1.tsv", sep = "\t", header = True, index = True)
target_raw_df.to_csv("data/tsv/aligned_raw_coords_sc2.tsv", sep = "\t", header = True, index = True)