#!/usr/bin/python

import os
import sys
from collections import defaultdict
import numpy as np
import utils_poca
import math
import json
import argparse
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import logging

log = utils_poca.Logger(__name__, level=logging.INFO)

# ========= Read the generated data ========= #

parser = argparse.ArgumentParser()
parser.add_argument('basename', type=str)
parser.add_argument('-r', '--recon', action='store_true')
parser.add_argument('-p', '--process', action='store_true')
parser.add_argument('-pp', '--post_process', action='store_true')
parser.add_argument('-v', '--voxel', type=float, default=10.0, help='Voxel size (default: 10.0)')
args = parser.parse_args()

basename = args.basename
do_reconstruction = args.recon
do_process = args.process
do_post_process = args.post_process
voxel_size = args.voxel


dataHandler = utils_poca.DataHandler([basename+".root"], basename, voxel_size=voxel_size)



# ========= Reconstruction ========= #
if do_reconstruction:
    dataHandler.load_data()
    dataHandler.reconstruction(test=False)  

# =========== Processing =========== #
if do_process:
    eventInfo_files = [f"{basename}_{dataHandler.voxel_size_str}_POCA_eventBasedVar_batch_*.npz"]
    eventInfo_files = utils_poca.get_list_make_list(eventInfo_files)
    if len(eventInfo_files) == 0:
       log.warning("No files found.")
    
    for eventInfo in eventInfo_files:
       dataHandler.read_eventInfo(eventInfo)
       dataHandler.update_S()
    dataHandler.update_lambda()
    
    density_storeTo = f"{basename}_{dataHandler.voxel_size_str}_POCA_scattering_density.npz"
    np.savez(density_storeTo, density=dataHandler.lambdaa)

# ========= Post-processing ========= #
if do_post_process:
    # Plot the signal voxels
    data = np.load(f"{basename}_{dataHandler.voxel_size_str}_POCA_scattering_density.npz")  # Load the .npz file
    dataHandler.lambdaa = data['density']
    dataHandler.plot_voxels()
    dataHandler.plot_voxels_3d()








