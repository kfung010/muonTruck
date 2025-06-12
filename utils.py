#!/usr/bin/python

import os
import sys
import uproot
import math
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from itertools import groupby, combinations
from operator import itemgetter
import multiprocessing
from functools import partial
from tqdm import tqdm

def Phi(vec):
    x, y, z = vec
    if x == 0 and y == 0:
        raise ValueError("Undefined phi angle for (0,0,z) vector")
    return math.atan2(y, x)

def toFilenameList(*names):  #Accepts any number of string arguments and returns them as a list.
    return list(names)

def checkFileExist(filenames):  #Check if the input files exist
    if not isinstance(filenames, list):
        filenames = [filenames]
    for fn in filenames:
        file_exists = os.path.exists(fn)
        if not file_exists:
            print(f"File '{fn}' does not exist.")
            sys.exit(1)

def replaceFile(filenames):  #Confirm to replace the existing output files
    if not isinstance(filenames, list):
        filenames = [filenames]
    for fn in filenames:
        outfile_exists = os.path.exists(fn)
        if outfile_exists:
            user_choice = input(f'\033[33mAre you sure to replace {fn}? [y/N] \033[0m').strip().lower()
            if user_choice != 'y':
                print("\033[31mOperation cancelled by user.\033[0m")
                sys.exit(1)

def returndf(filenames):

    if not isinstance(filenames, list):
        filenames = [filenames]
        
    hits_branches = "/eventNumber|muonMomentum|hitTime|hitPosition.|hitPixel.|hitMomentum./"
    #truck_branches = "/eventNumber|muonMomentum|scatTime|scatPosition.|scatMomentum./"    
    
    hits_dfs_grouped = []
    truck_dfs_grouped = []
    
    for fn in filenames:
    
        #print("Opening : ", fn )
        hits_tree = uproot.open(fn)["hits"]
        #truck_tree = uproot.open(fn)["truck"]
        #print("Done opening", fn)
        
        #print("  Grouping events ...")
        hits_b = hits_tree.keys(filter_name=hits_branches)
        #truck_b = truck_tree.keys(filter_name=truck_branches)
        
        hits_df = hits_tree.arrays(hits_b, library="pd")
        #truck_df = truck_tree.arrays(truck_b, library="pd")
        
        hits_df_grouped = hits_df.groupby("eventNumber")
        #truck_df_grouped = truck_df.groupby("eventNumber")
        #print("  Done grouping.")
        
        if len(filenames) == 1:
            return hits_df_grouped, []
            
        hits_dfs_grouped.append(hits_df_grouped)
        #truck_dfs_grouped.append(truck_df_grouped)
        
    return hits_dfs_grouped, truck_dfs_grouped

def getTotalEventNumber(df_grouped):
    return df_grouped.ngroups

def getBranchByEvent(df, eventNum, branch):
    if eventNum not in df.groups.keys():
        print(f"Event {eventNum} does not exist")
        return []
    return (df.get_group(eventNum))[branch].tolist()

def getBranchForAllEvents(df, branch):
    result = []
    total_events = getTotalEventNumber(df)
    '''
    for i, en in enumerate(df.groups.keys()):
        if (i % 100000 == 0):
            print(f"{i} events processed ... (Total: {total_events})")
        result.append(getBranchByEvent(df, en, branch))
    '''
    for i, (event_num, group) in enumerate(df):
        #if i % 100000 == 0:
        #    print(f"{i} events processed ... (Total: {total_events})")
        result.append(group[branch].tolist())
    return result

def process_chunk(chunk, branch):
    return [group[branch].tolist() for _, group in chunk]

def getBranchForAllEvents_parallel(df, branch, chunk_size=100000, num_processes=None):
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()
    
    groups = [(event_num, group) for event_num, group in df]
    total_events = len(groups)    
    chunks = [groups[i:i + chunk_size] for i in range(0, total_events, chunk_size)]
    
    process_chunk_with_branch = partial(process_chunk, branch=branch)
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        chunk_results = list(tqdm(
            pool.imap(process_chunk_with_branch, chunks),
            total=len(chunks),
            desc="Processing chunks"
        ))
    
    return [item for chunk in chunk_results for item in chunk]

def fit_line(x, y, z):  #Input the list of coordinates, fit a 3D line using SVD, and return the normalized direction

    points = np.vstack((x, y, z)).T

    centroid = np.mean(points, axis=0)  # Compute centroid
    _, _, Vt = np.linalg.svd(points - centroid)  # Singular Value Decomposition (SVD)
    direction = Vt[0]  # Extract principal direction
    return centroid, direction

def poca_reconstruction(hits_df, truck_df):

    results = {
        'hitXs_above': [],
        'hitYs_above': [],
        'hitZs_above': [],
        'hitXs_below': [],
        'hitYs_below': [],
        'hitZs_below': [],
        'numScat_list': [],
        'scatPositionX_truth_list': [],
        'scatPositionY_truth_list': [],
        'scatPositionZ_truth_list': [],
        'pocaXList': [],
        'pocaYList': [],
        'pocaZList': [],
        'pocaList': [],
        'above_centroid_list': [],
        'above_dirVec_list': [],
        'below_centroid_list': [],
        'below_dirVec_list': [],
        'angle_list': [],
        'proj_angle_list': [],
        'scat_dist' : []
    }
    
    total_events = getTotalEventNumber(hits_df)
    
    count = 0
    for (h_event, h_group) in hits_df:
    
        #if count % 100000 == 0:
        #    print(f"{count} events processed ... (Total: {total_events})")
        count += 1
        
        hitXs = h_group["hitPositionX_truth"].tolist()    
        hitYs = h_group["hitPositionY_truth"].tolist()
        hitZs = h_group["hitPositionZ_truth"].tolist()
              
        
        hitXs_above = [x for x, y in zip(hitXs, hitYs) if y > 0]
        hitXs_below = [x for x, y in zip(hitXs, hitYs) if y < 0]
        hitYs_above = [y for y in hitYs if y > 0]
        hitYs_below = [y for y in hitYs if y < 0]
        hitZs_above = [z for y, z in zip(hitYs, hitZs) if y > 0]
        hitZs_below = [z for y, z in zip(hitYs, hitZs) if y < 0]
        
        if len(hitYs_above) > 0 and len(hitYs_below) == 0:
            #print(f"Muon in event {h_event} is absorbed by the cargo. Skip this event.")
            continue

        if not (any(abs(x-y) > 100 for x, y in combinations(hitYs_above, 2)) and 
               any(abs(x-y) > 100 for x, y in combinations(hitYs_below, 2))):
            #print(f"Event {h_event} doesn't have enough RPC hits. Skipping.")
            continue

        above_centroid, above_dirVec = fit_line(hitXs_above, hitYs_above, hitZs_above)
        below_centroid, below_dirVec = fit_line(hitXs_below, hitYs_below, hitZs_below)
        
        w0 = above_centroid - below_centroid
        a = np.dot(above_dirVec, above_dirVec)
        b = np.dot(above_dirVec, below_dirVec)
        c = np.dot(below_dirVec, below_dirVec)
        d = np.dot(above_dirVec, w0)
        e = np.dot(below_dirVec, w0)
        denom = a * c - b * b
        
        if abs(denom) < 1e-11:
            #print(f"POCA reconstruction for event {h_event} is a straight line. Skipping.")
            continue
            
        s = (b * e - c * d) / denom
        t = (a * e - b * d) / denom
        
        poca1 = above_centroid + s * above_dirVec
        poca2 = below_centroid + t * below_dirVec
        poca = (poca1 + poca2) / 2
        
        scat_dist = math.sqrt(sum((a-b)**2 for a, b in zip(poca1, poca2)))
        dot_product = sum(a * b for a, b in zip(above_dirVec, below_dirVec))
        above_mom_mag = math.sqrt(sum(a**2 for a in above_dirVec))
        below_mom_mag = math.sqrt(sum(b**2 for b in below_dirVec))
        cos_theta_poca = dot_product / (above_mom_mag * below_mom_mag)
        angle_poca = math.acos(max(-1.0, min(1.0, cos_theta_poca))) * 1000  # mrad
        proj_angle = angle_poca * math.cos(math.atan2(below_dirVec[0], below_dirVec[2]))
        
        results['hitXs_above'].append(hitXs_above)
        results['hitYs_above'].append(hitYs_above)
        results['hitZs_above'].append(hitZs_above)
        results['hitXs_below'].append(hitXs_below)
        results['hitYs_below'].append(hitYs_below)
        results['hitZs_below'].append(hitZs_below)
        #results['numScat_list'].append(len(t_group["scatPositionX_truth"].tolist())-2)
        #results['scatPositionX_truth_list'].append(t_group["scatPositionX_truth"].tolist())
        #results['scatPositionY_truth_list'].append(t_group["scatPositionY_truth"].tolist())
        #results['scatPositionZ_truth_list'].append(t_group["scatPositionZ_truth"].tolist())
        results['pocaXList'].append(poca[0])
        results['pocaYList'].append(poca[1])   
        results['pocaZList'].append(poca[2])
        results['pocaList'].append(poca)
        results['above_centroid_list'].append(above_centroid)
        results['above_dirVec_list'].append(above_dirVec)
        results['below_centroid_list'].append(below_centroid)
        results['below_dirVec_list'].append(below_dirVec)
        results['angle_list'].append(angle_poca)
        results['proj_angle_list'].append(proj_angle)
        results['scat_dist'].append(scat_dist)

    return results
    
def create_box(center, size, color='gray', opacity=0.5):
    """Generate mesh3d vertices for a rectangular box."""
    x, y, z = center
    dx, dy, dz = size

    # Define the 8 corners
    x_vals = [x - dx / 2, x + dx / 2]
    y_vals = [y - dy / 2, y + dy / 2]
    z_vals = [z - dz / 2, z + dz / 2]

    # Vertices of the box
    vertices = np.array([
        [x_vals[0], y_vals[0], z_vals[0]],  # 0
        [x_vals[1], y_vals[0], z_vals[0]],  # 1
        [x_vals[1], y_vals[1], z_vals[0]],  # 2
        [x_vals[0], y_vals[1], z_vals[0]],  # 3
        [x_vals[0], y_vals[0], z_vals[1]],  # 4
        [x_vals[1], y_vals[0], z_vals[1]],  # 5
        [x_vals[1], y_vals[1], z_vals[1]],  # 6
        [x_vals[0], y_vals[1], z_vals[1]],  # 7
    ])
    
    # Define triangles for all 6 faces (each face has 2 triangles)
    triangles = [
        # Bottom face (z=z_vals[0])
        [0, 1, 2], [0, 2, 3],
        # Top face (z=z_vals[1])
        [4, 5, 6], [4, 6, 7],
        # Front face (y=y_vals[0])
        [0, 1, 5], [0, 5, 4],
        # Back face (y=y_vals[1])
        [3, 2, 6], [3, 6, 7],
        # Left face (x=x_vals[0])
        [0, 3, 7], [0, 7, 4],
        # Right face (x=x_vals[1])
        [1, 2, 6], [1, 6, 5],
    ]
    
    # Extract i, j, k indices for the triangles
    i = [t[0] for t in triangles]
    j = [t[1] for t in triangles]
    k = [t[2] for t in triangles]
    
    
    return go.Mesh3d(
        x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
        i=i, j=j, k=k,
        color=color, opacity=opacity
    )

def add_line_auto_range(fig, point, direction, color='purple', width=2):

    x_range = fig.layout.scene.xaxis.range
    y_range = fig.layout.scene.yaxis.range
    z_range = fig.layout.scene.zaxis.range
    
    t_values = []
    for axis, bounds in enumerate([x_range, y_range, z_range]):
        if direction[axis] != 0:
            t1 = (bounds[0] - point[axis]) / direction[axis]
            t2 = (bounds[1] - point[axis]) / direction[axis]
            t_values.extend([t1, t2])

    valid_t = [t for t in t_values if (
        (point[0] + direction[0]*t >= x_range[0]) and
        (point[0] + direction[0]*t <= x_range[1]) and
        (point[1] + direction[1]*t >= y_range[0]) and
        (point[1] + direction[1]*t <= y_range[1]) and
        (point[2] + direction[2]*t >= z_range[0]) and
        (point[2] + direction[2]*t <= z_range[1])
    )]
    
    if len(valid_t) >= 2:
        t_min = min(valid_t)
        t_max = max(valid_t)
    else:
        t_min, t_max = -150, 150 

    x_line = [point[0] + direction[0]*t_min, point[0] + direction[0]*t_max]
    y_line = [point[1] + direction[1]*t_min, point[1] + direction[1]*t_max]
    z_line = [point[2] + direction[2]*t_min, point[2] + direction[2]*t_max]

    fig.add_trace(go.Scatter3d(
        x=x_line,
        y=y_line,
        z=z_line,
        mode='lines',
        line=dict(color=color, width=width),
        name=f"Direction Line"
    )) 

def returnEventInfo(hits, truck, limit=-1):
    
    print("Reading RPC hit information ...")
    
    event_data = {}
    for event, x, y, z, mx, my, mz in zip(
        hits.eventNumber,
        hits.hitPositionX_truth,
        hits.hitPositionY_truth,
        hits.hitPositionZ_truth,
        hits.hitMomentumX_truth,
        hits.hitMomentumY_truth,
        hits.hitMomentumZ_truth
    ):      
        if event not in event_data:
            event_data[event] = {
                'hits_above': [], 'hits_below': [],
                'mom_above': [], 'mom_below': []
            }

        if y > 0:
            event_data[event]['hits_above'].append([x/10, z/10, y/10])
            event_data[event]['mom_above'].append([mx, mz, my])
        else:
            event_data[event]['hits_below'].append([x/10, z/10, y/10])
            event_data[event]['mom_below'].append([mx, mz, my])
    
    print("Done reading RPC hit information.")
    
    print("Reading truth scattering information ...")  
    scat_data = {}
    for event, x, y, z in zip(
        truck.eventNumber,
        truck.scatPositionX_truth,
        truck.scatPositionY_truth,
        truck.scatPositionZ_truth
    ):
        if event not in scat_data:
            scat_data[event] = []
        scat_data[event].append([x/10, z/10, y/10])    
    
    print("Done reading truth scattering information.")
    
    
    results = {
        'hitPosition_above_list': [],
        'hitPosition_below_list': [],
        'truth_scat_list': [],
        'above_centroid_list': [],
        'above_dirVec_list': [],
        'below_centroid_list': [],
        'below_dirVec_list': [],
        'poca_list': [],
        'angle_list': [],
        'scat_dist' : []
    }

    for event in sorted(event_data.keys()):
        if limit > 0 and len(results['poca_list']) >= limit:
            break
            
        data = event_data[event]
        hits_above = data['hits_above']
        hits_below = data['hits_below']
        
        if len(hits_above) > 0 and len(hits_below) == 0:
            print(f"Muon in event {event} is absorbed by the cargo. Skip this event.")
            continue

        z_above = [z for [_, _, z] in hits_above]
        z_below = [z for [_, _, z] in hits_below]
        
        if not (any(abs(x-y) > 3 for x, y in combinations(z_above, 2)) and 
               any(abs(x-y) > 3 for x, y in combinations(z_below, 2))):
            print(f"Event {event} doesn't have enough RPC hits. Skipping.")
            continue
            

        x_above, y_above, z_above = zip(*hits_above)
        x_below, y_below, z_below = zip(*hits_below)
        
        above_centroid, above_dirVec = fit_line(x_above, y_above, z_above)
        below_centroid, below_dirVec = fit_line(x_below, y_below, z_below)
        

        w0 = above_centroid - below_centroid
        a = np.dot(above_dirVec, above_dirVec)
        b = np.dot(above_dirVec, below_dirVec)
        c = np.dot(below_dirVec, below_dirVec)
        d = np.dot(above_dirVec, w0)
        e = np.dot(below_dirVec, w0)
        denom = a * c - b * b
        
        if abs(denom) < 1e-11:
            print(f"POCA reconstruction for event {event} is a straight line. Skipping.")
            continue
            
        s = (b * e - c * d) / denom
        t = (a * e - b * d) / denom
        
        poca1 = above_centroid + s * above_dirVec
        poca2 = below_centroid + t * below_dirVec
        poca = (poca1 + poca2) / 2
        
        scat_dist = math.sqrt(sum((a-b)**2 for a, b in zip(poca1, poca2)))
        dot_product = sum(a * b for a, b in zip(above_dirVec, below_dirVec))
        above_mom_mag = math.sqrt(sum(a**2 for a in above_dirVec))
        below_mom_mag = math.sqrt(sum(b**2 for b in below_dirVec))
        cos_theta_poca = dot_product / (above_mom_mag * below_mom_mag)
        angle_poca = math.acos(max(-1.0, min(1.0, cos_theta_poca))) * 1000  # mrad
        

        results['hitPosition_above_list'].append(hits_above)
        results['hitPosition_below_list'].append(hits_below)
        results['truth_scat_list'].append(scat_data.get(event, []))
        results['above_centroid_list'].append(above_centroid)
        results['above_dirVec_list'].append(above_dirVec)
        results['below_centroid_list'].append(below_centroid)
        results['below_dirVec_list'].append(below_dirVec)
        results['poca_list'].append(poca)
        results['angle_list'].append(angle_poca)
        results['scat_dist'].append(scat_dist)

    return (
        results['hitPosition_above_list'],
        results['hitPosition_below_list'],
        results['truth_scat_list'],
        np.array(results['above_centroid_list'], dtype=np.float64),
        np.array(results['above_dirVec_list'], dtype=np.float64),
        np.array(results['below_centroid_list'], dtype=np.float64),
        np.array(results['below_dirVec_list'], dtype=np.float64),
        np.array(results['poca_list'], dtype=np.float64),
        np.array(results['angle_list'], dtype=np.float64),
        np.array(results['scat_dist'], dtype=np.float64)
    )


















        