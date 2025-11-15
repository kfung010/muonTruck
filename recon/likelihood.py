#!/usr/bin/python

import os
import sys
from collections import defaultdict
import numpy as np
import utils_likelihood
import math
import json
import argparse
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import logging

log = utils_likelihood.Logger(__name__, level=logging.INFO)

# ========= Read the generated data and reconstruct the POCA ========= #

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

dataHandler = utils_likelihood.DataHandler([basename+".root"], basename, voxel_size=voxel_size)


# ========= Reconstruction ========= #
if do_reconstruction:
    dataHandler.load_data()
    dataHandler.reconstruction(test=False, useAve=False)

# =========== Processing =========== #
if do_process:
    eventInfo_files = [f"{basename}_{dataHandler.voxel_size_str}_likelihood_eventBasedVar_batch_*.npz"]
    eventInfo_files = utils_likelihood.get_list_make_list(eventInfo_files)
    if len(eventInfo_files) == 0:
        log.warning("No files found.")
    
    epochs = 50
    
    # Testing: plot the scattering density against epoch at some points
    n = 1
    epoch_arr = []
    lambda_vs_epoch = [[] for _ in range(n)]
    
    for epoch in range(epochs):
        log.info(f"Epoch {epoch}/{epochs}")
        total_files = len(eventInfo_files)
        for filenum, eventInfo in enumerate(eventInfo_files):
            percent = (filenum + 1) / total_files * 100
            filled = int(10 * (filenum + 1) // total_files)
            bar = '#' * filled + ' ' * (10 - filled)
            print(f"\rProcessing files: |{bar}| {percent:.1f}% ({filenum + 1} / {total_files})", end='', flush=True)
            dataHandler.read_eventInfo(eventInfo)

            # dataHandler.update_S_median()
            dataHandler.update_S()
            #dataHandler.update_S_angleOnly()
            #dataHandler.update_S_angleOnly_median()
            #dataHandler.update_S_dxOnly()
        
        #print()
        #print(dataHandler.times_being_hit[dataHandler.get_voxel_id(x=-0.1,y=0,z=0)])
        #print(dataHandler.times_being_hit[dataHandler.get_voxel_id(x=0.1,y=0,z=0)])

        
        # dataHandler.update_lambda_median()
        dataHandler.update_lambda()
        #dataHandler.update_lambda_angleOnly()
        #dataHandler.update_lambda_angleOnly_median()
        #dataHandler.update_lambda_dxOnly()
        # dataHandler.update_lambda_meanMedian()
        
        
        
        # Testing: print the lambda of the valid voxels
        valid = (dataHandler.lambdaa > 0.0000001) & (dataHandler.lambdaa < 1)
        valid_voxels = np.where(valid)[0]
        valid_voxels_lambda = dataHandler.lambdaa[valid]
        for idx, lam in zip(valid_voxels, valid_voxels_lambda):
            print(f"Index: {idx}, Lambda: {lam*1000000}")
        # Testing: plot the scattering density against epoch at some points
        epoch_arr.append(epoch+1)
        lambda_vs_epoch[0].append(dataHandler.lambdaa[dataHandler.get_voxel_id(x=-0.1,y=0,z=0)]*1000000)
        print()
    
    # Testing: plot the scattering density against epoch at some points
    # utils_likelihood.plot_line_chart(f"{basename}_lambda_vs_epoch_median.pdf", epoch_arr, lambda_vs_epoch, ["" for _ in range(n)], ["red" for _ in range(n)], "Epoch", "lambda")
    utils_likelihood.plot_line_chart(f"{basename}_lambda_vs_epoch_mean.pdf", epoch_arr, lambda_vs_epoch, ["" for _ in range(n)], ["red" for _ in range(n)], "Epoch", "lambda")
    
    
    dataHandler.lambdaa = dataHandler.lambdaa * 1000000
    density_storeTo = f"{basename}_{dataHandler.voxel_size_str}_likelihood_mean_scattering_density.npz"
    np.savez(density_storeTo, density=dataHandler.lambdaa)


# ========= Post-processing ========= #
if do_post_process:
    # Plot the signal voxels
    dataHandler.plot_voxels()
    dataHandler.plot_voxels_3d()






