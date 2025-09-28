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

basename = sys.argv[1]

dataHandler = utils_poca.DataHandler([basename+".root"], basename, voxel_size=10)
dataHandler.load_data()

# ========= Reconstruction ========= #
dataHandler.reconstruction(test=True)  # Comment out this line if it is already done

eventInfo_files = [f"{basename}_POCA_eventBasedVar_batch_*.npz"]
eventInfo_files = utils_poca.get_list_make_list(eventInfo_files)
if len(eventInfo_files) == 0:
    log.warning("No files found.")

for eventInfo in eventInfo_files:
    dataHandler.read_eventInfo(eventInfo)
    dataHandler.update_S()
dataHandler.update_lambda()

density_storeTo = f"{basename}_POCA_scattering_density.npz"
np.savez(density_storeTo, density=dataHandler.lambdaa)

# ========= Plot the signal voxels ========= #
dataHandler.plot_voxels()
dataHandler.plot_voxels_3d()








