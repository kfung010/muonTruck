#!/usr/bin/python

import os
import re
import sys
import math
import uproot
import json
import pprint
from collections import defaultdict
from functools import partial
import multiprocessing
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from scipy.stats import norm
from tqdm import tqdm
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from itertools import groupby, combinations
from operator import itemgetter
import logging
import glob

class Logger:
    VERBOSE_LEVEL_NUM = 15
    logging.addLevelName(VERBOSE_LEVEL_NUM, "VERBOSE")

    def __init__(self, name=__name__, level=logging.DEBUG):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)

        ch = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def verbose(self, message, *args, **kwargs):
        if self.logger.isEnabledFor(self.VERBOSE_LEVEL_NUM):
            self.logger._log(self.VERBOSE_LEVEL_NUM, message, args, **kwargs)

    def debug(self, message, *args, **kwargs):
        self.logger.debug(message, *args, **kwargs)

    def info(self, message, *args, **kwargs):
        self.logger.info(message, *args, **kwargs)

    def warning(self, message, *args, **kwargs):
        self.logger.warning(message, *args, **kwargs)

    def error(self, message, *args, **kwargs):
        self.logger.error(message, *args, **kwargs)

log = Logger(__name__, level=logging.INFO)

class CargoShape:
    def __init__(self, name, center, material):
        self.name = name
        self.center = center
        self.material = material
    
    def contains_point(self, x, y, z):
        # Return False by default, to be overwritten by derived class
        return False

class Box(CargoShape):
    def __init__(self, name, length, width, height, center, material):
        super().__init__(name, center, material)
        self.length = length
        self.width = width
        self.height = height
    def contains_point(self, x, y, z):
        cx, cy, cz = self.center
        return (cx - self.length/2 <= x <= cx + self.length/2 and
                cy - self.width/2  <= y <= cy + self.width/2  and
                cz - self.height/2 <= z <= cz + self.height/2)

class Ellipsoid(CargoShape):
    def __init__(self, name, rx, ry, rz, center, material):
        super().__init__(name, center, material)
        self.rx = rx
        self.ry = ry
        self.rz = rz
    def contains_point(self, x, y, z):
        cx, cy, cz = self.center
        dx = (x - cx) / self.rx
        dy = (y - cy) / self.ry
        dz = (z - cz) / self.rz
        return dx*dx + dy*dy + dz*dz <= 1

class Cylinder(CargoShape):
    def __init__(self, name, innerR, outerR, height, center, material):
        super().__init__(name, center, material)
        self.innerR = innerR
        self.outerR = outerR
        self.height = height
    def contains_point(self, x, y, z):
        cx, cy, cz = self.center
        dz = abs(z - cz)
        if dz > self.height / 2:
            return False
        dx = x - cx
        dy = y - cy
        r = (dx*dx + dy*dy)**0.5
        return self.innerR <= r <= self.outerR

def make_list(pattern):
    return glob.glob(pattern)

def get_list_make_list(input_data, input_dir = ""): 
    ''' Input a list of input data filenames (with wildcards), and output a full list of them '''
    if input_data == None:
        return []

    if not isinstance(input_data, list):
        input_data = [ input_data ]

    input_data_list = []

    for i, item in enumerate(input_data):
        if ( input_dir != "" ) and ( not item.startswith("/") ):
            input_data[i] = os.path.join(input_dir, item)
            item = input_data[i]
        input_data_list += make_list(item)

    return input_data_list


def generate_unique_filename(filename):
    if not os.path.exists(filename):
        return filename
    base, ext = os.path.splitext(filename)
    count = 1
    new_filename = f"{base}_{count}{ext}"
    while os.path.exists(new_filename):
        count += 1
        new_filename = f"{base}_{count}{ext}"
    return new_filename

def plot_histograms(outputPDF, datasets, labels, colors, lowLim, upLim, binnum, x_label, y_label, log_scale=False, figsize=(10,6), do_gauss_fit=False, fitLowlim=None, fitUplim=None):        
    ''' 
    Plot only one dataset:
    "datasets" is the list of data, e.g. [1.23, 4.56, ...]
    "labels" is a string (empty string means no label)
    "colors" is a string
    
    Plot multiple datasets:
    "datasets" is the list of the list of data, e.g. [[1.23, 4.56, ...], [1.34, 2.45, ...]]
    "labels" is the list of labels, e.g. ["x", "y"]
    "colors" is the list of colors, e.g. ["red", "green"]
    '''    
    
    if not isinstance(datasets[0], (list, np.ndarray, pd.Series)):
        datasets = [datasets]
    if isinstance(labels, str):
        labels = [labels] if labels else [None]
    if isinstance(colors, str):
        colors = [colors]
    
    bins = np.linspace(lowLim, upLim, binnum)
    
    outputPDF = generate_unique_filename(outputPDF)
    
    with PdfPages(outputPDF) as pdf:
        plt.figure(figsize=figsize)
        for data, label, color in zip(datasets, labels, colors):
            hist, bin_edges = np.histogram(data, bins=bins)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            if label:
                plt.step(bin_edges[:-1], hist, where='post', label=label, color=color)
            else:
                plt.step(bin_edges[:-1], hist, where='post', color=color)
            
            if do_gauss_fit:
                if fitLowlim is None:
                    fitLowlim = lowLim
                if fitUplim is None:
                    fitUplim = upLim
                try:
                    mask = (bin_centers >= fitLowlim) & (bin_centers <= fitUplim)
                    bin_centers_fit = bin_centers[mask]
                    hist_fit = hist[mask]
                    
                    
                    
                    
                    popt, _ = curve_fit(lambda x, mu, sigma, A: A * norm.pdf(x, mu, sigma),
                                    bin_centers_fit, hist_fit, p0=[np.mean(data), np.std(data), max(hist)])
                    mu, sigma, A = popt
                    x_fit_full = np.linspace(lowLim, upLim, 500)
                    y_fit_full = A * norm.pdf(x_fit_full, mu, sigma)
                    
                    mask_before = x_fit_full < fitLowlim
                    if np.any(mask_before):
                        plt.plot(x_fit_full[mask_before], y_fit_full[mask_before],
                                color=color, linestyle='--', linewidth=1.5)
                    
                    mask_in = (x_fit_full >= fitLowlim) & (x_fit_full <= fitUplim)
                    if np.any(mask_in):
                        fit_label = f'{label} Gaussian fit' if label else 'Gaussian fit'
                        plt.plot(x_fit_full[mask_in], y_fit_full[mask_in], 
                                color=color, linestyle='-', linewidth=1.5, label=fit_label)
                    
                    mask_after = x_fit_full > fitUplim
                    if np.any(mask_after):
                        plt.plot(x_fit_full[mask_after], y_fit_full[mask_after], 
                                color=color, linestyle='--', linewidth=1.5)

                    textstr = '\n'.join((
                        f'$\mu={mu:.5f}$',
                        f'$\sigma={sigma:.5f}$',
                    ))
                    xpos = mu
                    ypos = A * norm.pdf(mu, mu, sigma)
                    plt.text(xpos, ypos, textstr, fontsize=10,
                            color=color, verticalalignment='bottom', horizontalalignment='center',
                            bbox=dict(boxstyle='round,pad=0.3', edgecolor=color, facecolor='white', alpha=0.8))
                    
                except Exception as e:
                    print(f"Gaussian fit failed for dataset '{label}': {e}")
            
        plt.xlim(lowLim, upLim)
        plt.xlabel(x_label, fontsize=12)
        if log_scale:
            plt.yscale("log")
        plt.ylabel(y_label, fontsize=12)
        plt.legend()
        plt.grid(True, which='both', linestyle='--', alpha=0.5)
        pdf.savefig()
        plt.close()         
        
def plot_line_chart(outputPDF, x_values, y_values_list, labels, colors, x_label, y_label, linestyle='-', marker='o', figsize=(10, 6)):
    
    if not isinstance(y_values_list[0], (list, tuple)):
        y_values_list = [y_values_list]
    
    outputPDF = generate_unique_filename(outputPDF)
    
    with PdfPages(outputPDF) as pdf:
    
        plt.figure(figsize=figsize)
        
        for y_values, label, color in zip(y_values_list, labels, colors):
            
            plt.plot(x_values, y_values, label=label, color=color, linestyle=linestyle, marker=marker)
        
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.legend()
        
        plt.grid(True, linestyle='--', alpha=0.6)
        pdf.savefig()
        plt.close()     

def extract_largest_negative_smallest_positive(numbers):
    return max((n for n in numbers if n < 0), default=None), min((n for n in numbers if n > 0), default=None)

def fit_line(x, y, z):  
    '''Input the list of coordinates, fit a 3D line using SVD, and return the normalized direction'''
    points = np.vstack((x, y, z)).T
    centroid = np.mean(points, axis=0)  # Compute centroid
    _, _, Vt = np.linalg.svd(points - centroid)  # Singular Value Decomposition (SVD)
    direction = Vt[0]  # Extract principal direction
    return centroid, direction



# ========================================================================================

class DataHandler:

    def __init__(self, filenames, basename, configName=None, voxel_size=10):
    
        if not isinstance(filenames, list):
            filenames = [filenames]
        self.filenames = filenames
        self.basename = basename
        self.configName = configName

        self.p0 = 3000  # MeV

        self.recon_N_events = 0

        # Apparatus configuration
        self.pixelLength = None
        self.pixelNum1 = None
        self.pixelWidth = None
        self.pixelNum2 = None
        self.heights = []
        self.min_RPC_separation = None
        self.cargo_shapes = []
        self.x_min, self.x_max = None, None
        self.y_min, self.y_max = None, None
        self.z_min, self.z_max = None, None
        self.parse_config()

        # Voxelization configuration
        self.voxel_size = voxel_size
        self.x_min, self.x_max = None, None
        self.y_min, self.y_max = None, None
        self.z_min, self.z_max = None, None
        self.nx = None
        self.ny = None
        self.nz = None
        self.total_voxels = None
        self.voxel_size_str = str(voxel_size) + "cmVoxelSize"
        self.initialize_voxelization()

        # Simulation data
        self.hits_df_grouped = None

        # Per-event variables
        self.p_r2 = []

        # Per-voxel variables
        self.S_total = np.zeros(self.total_voxels)
        self.times_being_hit = np.zeros(self.total_voxels)
        self.lambdaa = np.full(self.total_voxels, 0.001)
        
        # =====

    def parse_config(self):
        
        if self.configName:
            try:
                with open(self.configName, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue
                        self._parse_line(line)
            except Exception as e:
                log.error(f"Error parsing {self.configName}: {e}")
        else:
            self.configName = self.basename+".mac"
            try:
                with open(self.configName, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue
                        self._parse_line(line)
            except Exception as e:
                log.error(f"Error parsing {self.configName}: {e}")
        
        if len(self.heights) > 0:
            min_separation = float('inf')  # Initialize with a large number
            for i in range(len(self.heights) - 1):
                separation = self.heights[i+1] - self.heights[i]
                if separation < min_separation:
                    min_separation = separation
            self.min_RPC_separation = min_separation

    def _parse_line(self, line):
    
        if line.startswith('/apparatus/pixelLength'):
            self.pixelLength = float(line.split()[1])
        elif line.startswith('/apparatus/pixelNum1'):
            self.pixelNum1 = float(line.split()[1])
        elif line.startswith('/apparatus/pixelWidth'):
            self.pixelWidth = float(line.split()[1])
        elif line.startswith('/apparatus/pixelNum2'):
            self.pixelNum2 = float(line.split()[1])
        elif line.startswith('/apparatus/addHeight'):
            heights = self.heights
            heights.append(float(line.split()[1]))
            self.heights = sorted(heights)
        elif line.startswith('/apparatus/addBox'):
            params = re.split(r'[,\s]+', line)[1:]
            box = Box(
                name=params[0],
                length=float(params[1]),
                width=float(params[2]),
                height=float(params[3]),
                center=(float(params[4]), float(params[5]), float(params[6])),
                material=params[7]
            )
            self.cargo_shapes.append(box)

        elif line.startswith('/apparatus/addEllipsoid'):
            params = re.split(r'[,\s]+', line)[1:]
            ellipsoid = Ellipsoid(
                name=params[0],
                rx=float(params[1]),
                ry=float(params[2]),
                rz=float(params[3]),
                center=(float(params[4]), float(params[5]), float(params[6])),
                material=params[7]
            )
            self.cargo_shapes.append(ellipsoid)

        elif line.startswith('/apparatus/addCylinder'):
            params = re.split(r'[,\s]+', line)[1:]
            cylinder = Cylinder(
                name=params[0],
                innerR=float(params[1]),
                outerR=float(params[2]),
                height=float(params[3]),
                center=(float(params[4]), float(params[5]), float(params[6])),
                material=params[7]
            )
            self.cargo_shapes.append(cylinder)

    def initialize_voxelization(self):

        log.info(f"Using {self.voxel_size} cm voxelization.")

        total_length = self.pixelLength * self.pixelNum1
        total_width = self.pixelWidth * self.pixelNum2
        
        if (self.x_min is None) or (self.y_min is None) or (self.z_min is None):
            self.x_min, self.x_max = -total_length / 2 , total_length / 2
            self.y_min, self.y_max = -total_width / 2 , total_width / 2
            self.z_min, self.z_max = extract_largest_negative_smallest_positive(self.heights)
        
        self.nx = int((self.x_max - self.x_min) / self.voxel_size)
        self.ny = int((self.y_max - self.y_min) / self.voxel_size)
        self.nz = int((self.z_max - self.z_min) / self.voxel_size)
        
        self.total_voxels = self.nx * self.ny * self.nz

    def load_data(self):

        log.info("Loading data ...")
    
        hits_branches = "/eventNumber|muonMomentum|hitTime|hitPosition.|hitPixel.|hitMomentum./"
        
        all_hits_dfs = []
        
        for fn in self.filenames:
            hits_tree = uproot.open(fn)["hits"]
            hits_b = hits_tree.keys(filter_name=hits_branches)
            hits_df = hits_tree.arrays(hits_b, library="pd")
            all_hits_dfs.append(hits_df)
        
        combined_hits_df = pd.concat(all_hits_dfs, ignore_index=True)    
        self.hits_df_grouped = combined_hits_df.groupby("eventNumber")

        log.info("Loaded data.")

    def read_eventInfo(self, data_filename):
        log.debug(f"Loading {data_filename} ...")
        data = np.load(data_filename, allow_pickle=True)
        self.p_r2 = data['p_r2']
        self.D_x = data['D_x']
        self.D_y = data['D_y']
        self.passed_voxels = data['passed_voxels']
        self.POCA_voxel = data["POCA_voxel"]
        log.debug(f"Loaded {data_filename}.")
    
    def reconstruction(self, test=False, batch_size=1000, suffix=""):
        
        log.info("Reconstruction starts.")

        if self.hits_df_grouped is None:
            self.load_data()

        total_events = len(self.hits_df_grouped)

        batch_count = 0
        batch_p_r2 = []
        batch_D_x = []
        batch_D_y = []
        batch_S_x = []
        batch_passed_voxels = []
        batch_POCA_voxel = []

        test_num = 1000

        for count, (h_event, h_group) in enumerate(self.hits_df_grouped):
            if test and count > test_num:
                break

            percent = (count + 1) / total_events * 100
            filled = int(10 * (count + 1) // total_events)
            bar = '#' * filled + ' ' * (10 - filled)
            print(f"\rProcessing events: |{bar}| {percent:.1f}% ({count + 1} / {total_events}, {self.recon_N_events} reconstructed)", end='', flush=True)


            # ***********************************************************************************
            # Fit two straight lines for the RPC hits (above and below the region of interest)
            # ***********************************************************************************
            hits = np.array([[x/10, y/10, z/10] for x, y, z in 
                                zip(h_group["hitPositionX_truth"].tolist(),
                                    h_group["hitPositionY_truth"].tolist(),
                                    h_group["hitPositionZ_truth"].tolist())])

            hits_above = hits[hits[:, 2] > 0]
            hits_below = hits[hits[:, 2] < 0]

            # Skip events that are absorbed by the cargo
            if len(hits_above) > 0 and len(hits_below) == 0:
                continue
            
            # Skip events that do not have engough RPC hits
            def _has_large_gap(z_list):
                return any(abs(a - b) > self.min_RPC_separation/2 for a, b in combinations(z_list, 2))
            if not (_has_large_gap([z for _, _, z in hits_above]) and 
                    _has_large_gap([z for _, _, z in hits_below])):
                continue

            above_centroid, above_dirVec = fit_line(*zip(*hits_above))
            below_centroid, below_dirVec = fit_line(*zip(*hits_below))
            
            # ***********************************************************************************
            # POCA (point-of-closest-approach) reconstruction
            # ***********************************************************************************
            w0 = above_centroid - below_centroid
            a = np.dot(above_dirVec, above_dirVec)
            b = np.dot(above_dirVec, below_dirVec)
            c = np.dot(below_dirVec, below_dirVec)
            d = np.dot(above_dirVec, w0)
            e = np.dot(below_dirVec, w0)
            denom = a * c - b * b
            
            # Skip events where the POCA reconstruction is essentially a straight line
            if abs(denom) < 1e-11:
                continue
  
            s = (b * e - c * d) / denom
            t = (a * e - b * d) / denom
            
            poca1 = above_centroid + s * above_dirVec
            poca2 = below_centroid + t * below_dirVec
            poca = (poca1 + poca2) / 2

            # Skip events where the POCA is not in the region of interest
            if (poca[0] < self.x_min or poca[0] > self.x_max or 
                poca[1] < self.y_min or poca[1] > self.y_max or 
                poca[2] < self.z_min or poca[2] > self.z_max):
                continue

            # Relative momentum squared
            momentum_vector = h_group[["hitMomentumX_truth", "hitMomentumY_truth", "hitMomentumZ_truth"]].iloc[0].values
            momentum = np.linalg.norm(momentum_vector)
            p_r2 = (self.p0 / (momentum*1000))**2

            # Entry and exit point of muon            
            point_0 = self.intersections(above_centroid, above_dirVec, self.z_max)
            point_p = self.intersections(above_centroid, above_dirVec, self.z_min, requireInRange=False)
            point_1 = self.intersections(below_centroid, below_dirVec, self.z_min)        

            if (point_0 is None) or (point_p is None) or (point_1 is None):
                continue

            # Incident angles and projected scattering angles
            theta_x_0, theta_y_0, delta_theta_x, delta_theta_y = DataHandler.angles(above_dirVec, below_dirVec)
            
            # Displacements
            Lxy = math.sqrt(1 + math.tan(theta_x_0) ** 2 + math.tan(theta_y_0) ** 2)
            delta_x = (point_1[0]-point_p[0]) * math.cos(theta_x_0) * Lxy * math.cos(delta_theta_x+theta_x_0) / math.cos(delta_theta_x)
            delta_y = (point_1[1]-point_p[1]) * math.cos(theta_y_0) * Lxy * math.cos(delta_theta_y+theta_y_0) / math.cos(delta_theta_y)

            # Data vector
            D_x = delta_theta_x*1000
            D_y = delta_theta_y*1000

            # Voxels passed by the POCA-reconstructed-trajectory (entry-POCA-exit)
            passed_voxels_before_POCA = self.get_voxels_along_line(point_0, poca)
            passed_voxels_after_POCA = self.get_voxels_along_line(poca, point_1)
            passed_voxels = np.concatenate((passed_voxels_before_POCA[:-1], passed_voxels_after_POCA))
            POCA_voxel = self.get_voxel_id(x=poca[0],y=poca[1],z=poca[2])

            self.recon_N_events += 1

            batch_p_r2.append(p_r2)
            batch_D_x.append(D_x)
            batch_D_y.append(D_y)
            batch_passed_voxels.append(passed_voxels)
            batch_POCA_voxel.append(POCA_voxel)

            stop = (count + 1) % batch_size == 0 or (count + 1) == total_events
            if test:
                stop = (count + 1) % batch_size == 0 or count == test_num
            if stop:
                batch_filename = f"{self.basename}_{self.voxel_size_str}_POCA_eventBasedVar_batch_{batch_count}.npz"
                np.savez(
                    batch_filename,
                    p_r2=np.array(batch_p_r2),
                    D_x=np.array(batch_D_x),
                    D_y=np.array(batch_D_y),
                    passed_voxels=np.array(batch_passed_voxels, dtype=object),
                    POCA_voxel=np.array(batch_POCA_voxel)
                )

                batch_count += 1

                batch_p_r2 = []
                batch_D_x = []
                batch_D_y = []
                batch_passed_voxels = []
                batch_POCA_voxel = []

        print()
        log.info("Reconstruction ends.")

    def update_S(self):
        
        for count, (p_r2, D_x, D_y, passed_voxels, POCA_voxel) in enumerate(zip(self.p_r2, self.D_x, self.D_y, self.passed_voxels, self.POCA_voxel)):
            
            self.times_being_hit[passed_voxels] += 1
            self.S_total[POCA_voxel] += (D_x*D_x + D_y*D_y) / 2 

    def update_lambda(self):
        valid = self.times_being_hit != 0
        self.lambdaa[valid] = (self.S_total[valid] / self.times_being_hit[valid]) / self.voxel_size / 5.1526


    # ===============================================================================

    def get_voxel_id(self, x=None, y=None, z=None, indices=None):
        if indices is not None:
            ix, iy, iz = indices
        else:
            ix = int((x - self.x_min) // self.voxel_size)
            iy = int((y - self.y_min) // self.voxel_size)
            iz = int((z - self.z_min) // self.voxel_size)
            ix = max(0, min(ix, self.nx - 1))
            iy = max(0, min(iy, self.ny - 1))
            iz = max(0, min(iz, self.nz - 1))
        return (iz * self.ny + iy) * self.nx + ix
    
    def get_voxel_indices(self, voxel_id=None, x=None, y=None, z=None):
        if voxel_id is not None:
            ix = voxel_id % self.nx
            iy = (voxel_id // self.nx) % self.ny
            iz = voxel_id // (self.nx * self.ny)
        else:
            ix = int((x - self.x_min) // self.voxel_size)
            iy = int((y - self.y_min) // self.voxel_size)
            iz = int((z - self.z_min) // self.voxel_size)
            ix = max(0, min(ix, self.nx - 1))
            iy = max(0, min(iy, self.ny - 1))
            iz = max(0, min(iz, self.nz - 1))
        return np.array([ix, iy, iz])
    
    def get_voxel_coor(self, voxel_id=None, indices=None):
        if voxel_id is not None:
            indices = self.get_voxel_indices(voxel_id=voxel_id)
        if indices is not None:
            ix, iy, iz = indices
            x = self.x_min + (ix + 0.5) * self.voxel_size
            y = self.y_min + (iy + 0.5) * self.voxel_size
            z = self.z_min + (iz + 0.5) * self.voxel_size
            return np.array([x, y, z])
        return None

    def intersections(self, point, dirVec, height, requireInRange=True):
        x0, y0, z0 = point
        vx, vy, vz = dirVec
        if vz != 0:
            t = (height - z0) / vz
            x = x0 + t * vx
            y = y0 + t * vy
            z = height
            if requireInRange:
                if x >= self.x_min and x <= self.x_max and y > self.y_min and y <= self.y_max:
                    return np.array([x, y, z])
                else:
                    return None
            else:
                return np.array([x, y, z])
        else:
            return None

    @classmethod
    def angles(cls, above_dirVec, below_dirVec):
        entry_angle_xz = math.atan(above_dirVec[0]/abs(above_dirVec[2]))
        entry_angle_yz = math.atan(above_dirVec[1]/abs(above_dirVec[2]))
        exit_angle_xz = math.atan(below_dirVec[0]/abs(below_dirVec[2]))
        exit_angle_yz = math.atan(below_dirVec[1]/abs(below_dirVec[2]))
        d_theta_x = exit_angle_xz - entry_angle_xz
        d_theta_y = exit_angle_yz - entry_angle_yz
        return entry_angle_xz, entry_angle_yz, d_theta_x, d_theta_y
    
    def passed_voxel_min_max(self, passed_voxels):
        ix = passed_voxels % self.nx
        iy = (passed_voxels // self.nx) % self.ny
        iz = passed_voxels // (self.nx * self.ny)
        x_coords = self.x_min + (ix + 0.5) * self.voxel_size
        y_coords = self.y_min + (iy + 0.5) * self.voxel_size
        z_coords = self.z_min + (iz + 0.5) * self.voxel_size
        voxel_center = np.stack([x_coords, y_coords, z_coords], axis=1)
        voxel_min = voxel_center - self.voxel_size / 2.0
        voxel_max = voxel_center + self.voxel_size / 2.0
        return voxel_min, voxel_max

    def liang_barsky(self, p1, p2, voxel_mins, voxel_maxs):

        # Direction vector from p1 to p2
        d = p2 - p1
        voxel_mins = np.array(voxel_mins)
        voxel_maxs = np.array(voxel_maxs)

        # Initialize parameter range to [0,1]
        tmin = np.zeros(len(voxel_mins))
        tmax = np.ones(len(voxel_mins))

        # Process each dimension (x, y, z) separately
        for i in range(3):
            di = d[i]   # Direction component
            p1i = p1[i]   # Start point component
            min_i = voxel_mins[:, i]
            max_i = voxel_maxs[:, i]

            # Handel case when line segment is parallel to voxel face
            parallel_mask = (di == 0)
            outside_parallel = parallel_mask & ((p1i < min_i) | (p1i > max_i))
            if np.any(outside_parallel):
                tmin[outside_parallel] = np.nan
                tmax[outside_parallel] = np.nan
            
            # Non-parallel case
            non_parallel_mask = ~parallel_mask
            t1 = np.full(len(voxel_mins), -np.inf)  # Parameter for first boundary
            t2 = np.full(len(voxel_mins), np.inf)   # Parameter for second boundary

            # Calculate intersection parameters with voxel boundaries
            t1[non_parallel_mask] = (min_i[non_parallel_mask] - p1i) / di
            t2[non_parallel_mask] = (max_i[non_parallel_mask] - p1i) / di

            # Determine which boundary is entry (t_low) and which is exit (t_high)
            t_low = np.minimum(t1, t2)
            t_high = np.maximum(t1, t2)

            # Update global parameter range
            tmin = np.maximum(tmin, t_low)
            tmax = np.minimum(tmax, t_high)

        epsilon = 1e-4
        nearly_equal = np.abs(tmin - tmax) < epsilon
        avg_value = (tmin + tmax) / 2
        tmin = np.where(nearly_equal, avg_value, tmin)
        tmax = np.where(nearly_equal, avg_value, tmax)
        valid = (~np.isnan(tmin)) & (~np.isnan(tmax)) & (tmin <= tmax) & (tmax >= 0) & (tmin <= 1)
        tmin_clipped = np.where(valid, np.clip(tmin, 0, 1), np.nan)
        tmax_clipped = np.where(valid, np.clip(tmax, 0, 1), np.nan)

        return tmin_clipped, tmax_clipped

    def get_voxels_along_line(self, p0, p1):

        start_voxel = self.get_voxel_indices(x=p0[0], y=p0[1], z=p0[2])
        end_voxel = self.get_voxel_indices(x=p1[0], y=p1[1], z=p1[2])
        
        if np.array_equal(start_voxel, end_voxel):
            return np.array([self.get_voxel_id(indices=start_voxel)])
        
        # Calculate the t parameter at the boundary of the next voxel for each direction
        direction = p1 - p0
        step = np.sign(direction).astype(int)
        t_max = np.zeros(3)
        for i in range(3):
            if direction[i] != 0:
                axis_min = getattr(self, f"{['x', 'y', 'z'][i]}_min")
                current_voxel_min = axis_min + start_voxel[i] * self.voxel_size
                current_voxel_max = current_voxel_min + self.voxel_size
                if step[i] > 0:
                    next_boundary = current_voxel_max
                else:
                    next_boundary = current_voxel_min
                t_max[i] = (next_boundary - p0[i]) / direction[i]
            else:
                t_max[i] = float('inf')
        
        # Calculate the change of t for a unit step along the line segment
        t_delta = np.zeros(3)
        for i in range(3):
            if direction[i] != 0:
                t_delta[i] = self.voxel_size / abs(direction[i])
            else:
                t_delta[i] = float('inf')
        
        # Iteratively find the axes with smallest t_max
        current_voxel = start_voxel.copy()
        voxel_ids = [self.get_voxel_id(indices=current_voxel)]
        max_steps = 2 * (self.nx + self.ny + self.nz)
        steps = 0
        tolerance = 1e-6
        while steps < max_steps:

            steps += 1
            min_t = np.min(t_max)
            axes_with_min_t = []
            for i in range(3):
                if abs(t_max[i] - min_t) < tolerance:
                    axes_with_min_t.append(i)
            
            beyond_end = False
            for i in range(3):
                if step[i] > 0 and current_voxel[i] > end_voxel[i]:
                    beyond_end = True
                    break
                elif step[i] < 0 and current_voxel[i] < end_voxel[i]:
                    beyond_end = True
                    break
            if beyond_end:
                break

            if len(axes_with_min_t) > 1:   # Multiple smallest t_max
                moved_axes = []
                for axis in axes_with_min_t:
                    current_voxel[axis] += step[axis]
                    moved_axes.append(axis)
                    if (current_voxel[0] < 0 or current_voxel[0] >= self.nx or
                        current_voxel[1] < 0 or current_voxel[1] >= self.ny or
                        current_voxel[2] < 0 or current_voxel[2] >= self.nz):
                        break
                    t_max[axis] += t_delta[axis]
                    voxel_id = self.get_voxel_id(indices=current_voxel)
                    voxel_ids.append(voxel_id)
                    if np.array_equal(current_voxel, end_voxel):
                        break
                if np.array_equal(current_voxel, end_voxel):
                    break
            else:
                axis = axes_with_min_t[0]
                current_voxel[axis] += step[axis]
                if (current_voxel[0] < 0 or current_voxel[0] >= self.nx or
                    current_voxel[1] < 0 or current_voxel[1] >= self.ny or
                    current_voxel[2] < 0 or current_voxel[2] >= self.nz):
                    break
                t_max[axis] += t_delta[axis]
                voxel_id = self.get_voxel_id(indices=current_voxel)
                voxel_ids.append(voxel_id)
                if np.array_equal(current_voxel, end_voxel):
                    break
        
        return np.array(voxel_ids)

    def generate_shapes_for_view(self, view):
        view_shapes = []
        for shape in self.cargo_shapes:
            if isinstance(shape, Box):
                cx, cy, cz = shape.center
                if view == 'xy':
                    view_shapes.append({
                        'type': 'rect', 'center': (cx, cy),
                        'length': shape.length,'width': shape.width,
                        'color': 'red'
                    })
                elif view == 'xz':
                    view_shapes.append({
                        'type': 'rect', 'center': (cx, cz),
                        'length': shape.length, 'width': shape.height,
                        'color': 'red'
                    })
                elif view == 'yz':
                    view_shapes.append({
                        'type': 'rect', 'center': (cy, cz),
                        'length': shape.width, 'width': shape.height,
                        'color': 'red'
                    })
    
            elif isinstance(shape, Ellipsoid):
                cx, cy, cz = shape.center
                if view == 'xy':
                    view_shapes.append({
                        'type': 'ellipse', 'center': (cx, cy),
                        'length': 2 * shape.rx, 'width': 2 * shape.ry,
                        'color': 'red'
                    })
                elif view == 'xz':
                    view_shapes.append({
                        'type': 'ellipse', 'center': (cx, cz),
                        'length': 2 * shape.rx, 'width': 2 * shape.rz,
                        'color': 'red'
                    })
                elif view == 'yz':
                    view_shapes.append({
                        'type': 'ellipse', 'center': (cy, cz),
                        'length': 2 * shape.ry, 'width': 2 * shape.rz,
                        'color': 'red'
                    })
    
            elif isinstance(shape, Cylinder):
                cx, cy, cz = shape.center
                if view == 'xy':
                    view_shapes.append({
                        'type': 'circle', 'center': (cx, cy),
                        'radius': shape.outerR,
                        'color': 'red'
                    })
                elif view == 'xz':
                    view_shapes.append({
                        'type': 'rect', 'center': (cx, cz),
                        'length': 2 * shape.outerR, 'width': shape.height,
                        'color': 'red'
                    })
                elif view == 'yz':
                    view_shapes.append({
                        'type': 'rect', 'center': (cy, cz),
                        'length': 2 * shape.outerR, 'width': shape.height,
                        'color': 'red'
                    })
    
        return view_shapes


    def plot_voxels(self, classified=True, lowLim=0.1, upLim=1000000, cmap='rainbow', alpha=0.7):
        
        valid = (self.lambdaa > lowLim) & (self.lambdaa < upLim)
        valid_voxels = np.where(valid)[0]
        valid_voxels_lambda = self.lambdaa[valid]

        if len(valid_voxels) == 0:
            log.warning("No valid voxels to plot.")
            return
        
        if not classified:
            min_val, max_val = min(valid_voxels_lambda), max(valid_voxels_lambda)
            norm = plt.Normalize(vmin=min_val, vmax=max_val)
            cmap = plt.get_cmap(cmap)
        
        rpc_length = self.x_max - self.x_min
        rpc_width = self.y_max - self.y_min
        z_min, z_max = self.z_min, self.z_max
        
        views = ['xy', 'xz', 'yz']
        
        outputPDF = generate_unique_filename(f"valid_voxels_projected_poca_{self.basename}_{self.voxel_size_str}.pdf")
        with PdfPages(outputPDF) as pdf:
            for view in views:
                fig, ax = plt.subplots(figsize=(6, 6))
                if view == 'xy':
                    ax.set_xlabel('x [cm]')
                    ax.set_ylabel('y [cm]')
                    ax.set_xlim(self.x_min, self.x_max)
                    ax.set_ylim(self.y_min, self.y_max)
                elif view == 'xz':
                    ax.set_xlabel('x [cm]')
                    ax.set_ylabel('z [cm]')
                    ax.set_xlim(self.x_min, self.x_max)
                    ax.set_ylim(z_min, z_max)
                else:  # yz
                    ax.set_xlabel('y [cm]')
                    ax.set_ylabel('z [cm]')
                    ax.set_xlim(self.y_min, self.y_max)
                    ax.set_ylim(z_min, z_max)
                ax.set_aspect('equal')
                ax.grid(True, linestyle='--', alpha=0.5)
                
                for vid, lambdaa in zip(valid_voxels, valid_voxels_lambda):
                    
                    if classified:
                        if 0.1 < lambdaa < 5:
                            color = 'pink'
                        elif 5 < lambdaa < 30:
                            color = 'green'
                        elif lambdaa > 30:
                            color = 'red'
                        else:
                            color = 'grey'
                    else:
                        color = cmap(norm(lambdaa))
                    
                    cx, cy, cz = self.get_voxel_coor(voxel_id=vid)
                    voxel_size = self.voxel_size
                    if view == 'xy':
                        center = (cx, cy)
                        length, width = voxel_size, voxel_size
                    elif view == 'xz':
                        center = (cx, cz)
                        length, width = voxel_size, voxel_size
                    else:  # yz
                        center = (cy, cz)
                        length, width = voxel_size, voxel_size
                    
                    lower_left_x = center[0] - length / 2
                    lower_left_y = center[1] - width / 2
                    rect = patches.Rectangle((lower_left_x, lower_left_y), length, width,
                                             linewidth=0, edgecolor=None,
                                             facecolor=color, alpha=alpha)
                    ax.add_patch(rect)
                
                shapes_for_view = self.generate_shapes_for_view(view)
                for shape in shapes_for_view:
                    if shape["type"] == "rect":
                        cx, cy = shape["center"]
                        length = shape["length"]
                        width = shape["width"]
                        lower_left_x = cx - length / 2
                        lower_left_y = cy - width / 2
                        rect = patches.Rectangle(
                            (lower_left_x, lower_left_y), length, width,
                            linewidth=2, edgecolor=shape["color"], facecolor='none'
                        )
                        ax.add_patch(rect)
                    elif shape["type"] == "ellipse":
                        cx, cy = shape["center"]
                        length = shape["length"]
                        width = shape["width"]
                        ellipse = patches.Ellipse(
                            (cx, cy), length, width,
                            linewidth=2, edgecolor=shape["color"], facecolor='none'
                        )
                        ax.add_patch(ellipse)
                    elif shape["type"] == "circle":
                        cx, cy = shape["center"]
                        radius = shape["radius"]
                        circle = patches.Circle(
                            (cx, cy), radius,
                            linewidth=2, edgecolor=shape["color"], facecolor='none'
                        )
                        ax.add_patch(circle)
                
                if not classified:    
                    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                    sm.set_array([])
                    cbar = fig.colorbar(sm, ax=ax)
                    cbar.set_label("Scattering density [1/cm]")
                
                pdf.savefig(fig)
                plt.close(fig)


    def plot_voxels_3d(self, classified=True, lowLim=0.1, upLim=1000000, cmap='rainbow', alpha=0.7):
        
        valid = (self.lambdaa > lowLim) & (self.lambdaa < upLim)
        valid_voxels = np.where(valid)[0]
        valid_voxels_lambda = self.lambdaa[valid]

        if len(valid_voxels) == 0:
            log.warning("No valid voxels to plot.")
            return
        
        x_coords = []
        y_coords = []
        z_coords = []
        colors = []
        
        if not classified:
            min_val, max_val = min(valid_voxels_lambda), max(valid_voxels_lambda)
            norm = plt.Normalize(vmin=min_val, vmax=max_val)
            cmap = plt.get_cmap(cmap)
        
        for vid, lambdaa in zip(valid_voxels, valid_voxels_lambda):
            x, y, z = self.get_voxel_coor(voxel_id=vid)
            x_coords.append(x)
            y_coords.append(y)
            z_coords.append(z)
            
            if classified:
                if 0.1 < lambdaa < 5:
                    colors.append('pink')
                elif 5 < lambdaa < 30:
                    colors.append('green')
                elif lambdaa > 30:
                    colors.append('red')
                else:
                    colors.append('grey')
            else:
                colors.append(lambdaa)
        
        fig = go.Figure()
        print()
        if classified:
            fig.add_trace(go.Scatter3d(
                x=x_coords,
                y=y_coords,
                z=z_coords,
                mode='markers',
                marker=dict(
                    size=3, 
                    color=colors,
                    opacity=alpha,
                ),
            ))
        else:
            fig.add_trace(go.Scatter3d(
                x=x_coords,
                y=y_coords,
                z=z_coords,
                mode='markers',
                marker=dict(
                    size=3, 
                    color=colors,
                    colorscale=cmap, 
                    opacity=alpha, 
                    colorbar=dict(title="Scattering density [1/cm]")
                ),
            ))
        
        low_bound = min(self.x_min, self.y_min, self.z_min)
        up_bound = max(self.x_max, self.y_max, self.z_max)

        fig.update_layout(
            title="",
            scene=dict(
                xaxis=dict(title="x (cm)", range=[low_bound, up_bound]),
                yaxis=dict(title="y (cm)", range=[low_bound, up_bound]),
                zaxis=dict(title="z (cm)", range=[low_bound, up_bound]),
                aspectmode="cube"
            ),
            margin=dict(l=0, r=0, t=40, b=0),
        )
        
        outputHTML = generate_unique_filename(f"valid_voxels_poca_{self.basename}_{self.voxel_size_str}.html")
        fig.write_html(outputHTML)

    
    def is_voxel_truth_signal(self, voxel_id):
        x, y, z = self.get_voxel_coor(voxel_id=voxel_id)
        isSignal = False
        for shape in self.cargo_shapes:
            if shape.contains_point(x, y, z):
                isSignal = True
                break
        return isSignal

    def evaluate(self):
        thresholds = np.arange(0.1, 2.1, 0.1)
        accuracies = []
        precisions = []
        recalls = []
        f1_scores = []
        false_rates = []

        truth = np.zeros(self.total_voxels)
        for voxel_id in range(self.total_voxels):
            if self.is_voxel_truth_signal(voxel_id):
                truth[voxel_id] = 1

        for threshold in thresholds:
            recon = np.where(self.lambdaa < threshold, 0, 1)

            N1 = np.sum((truth == 1) & (recon == 1))
            N2 = np.sum((truth == 0) & (recon == 1))
            N3 = np.sum((truth == 1) & (recon == 0))
            N4 = np.sum((truth == 0) & (recon == 0))

            #     Accuracy : Fraction of correctly identified voxels
            #     Precision : Fraction of predicted signal voxels that are actually signal
            #     Recall : Fraction of actual signal voxels that are predicted signal
            #     F1 Score : Harmonic mean of precision and recall
            #     False rate : Fraction of actual background voxels wrongly predicted as signal
            accuracy = (N1 + N4) / (N1 + N2 + N3 + N4) if (N1 + N2 + N3 + N4) > 0 else 0
            precision = N1 / (N1 + N2) if (N1 + N2) > 0 else 0
            recall = N1 / (N1 + N3) if (N1 + N3) > 0 else 0
            f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
            false_rate = N2 / (N4 + N2) if (N4 + N2) > 0 else 0

            accuracies.append(accuracy)
            precisions.append(precision)
            recalls.append(recall)
            f1_scores.append(f1_score)
            false_rates.append(false_rate)

        outputPDF = generate_unique_filename(f"metric_poca_{self.basename}_{self.voxel_size_str}.pdf")
        with PdfPages(outputPDF) as pdf:
            metrics = {
                'Accuracy': accuracies,
                'Precision': precisions,
                'Recall': recalls,
                'F1 Score': f1_scores,
                'False Rate': false_rates
            }

            for metric_name, values in metrics.items():
                fig, ax = plt.subplots(figsize=(8, 6))
                ax.plot(thresholds, values, marker='o')
                ax.set_xlabel('Scattering density cut')
                ax.set_ylabel(metric_name)
                ax.set_xlim(0.1, 2.0)
                if metric_name == "Accuracy":
                    ax.set_ylim(0.9, 1)
                elif metric_name == "False Rate":
                    ax.set_ylim(0.0, 0.1)
                else:
                    ax.set_ylim(0., 1.)
                
                ax.grid(True)

                pdf.savefig(fig)
                plt.close(fig)