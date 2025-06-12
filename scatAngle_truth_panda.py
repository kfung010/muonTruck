#!/usr/bin/python

import os
import sys
import uproot
import pandas as pd
import utils
import numpy as np
import plotly.graph_objects as go
from itertools import groupby
from operator import itemgetter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from matplotlib import gridspec

# ==================================================================


# ==================================================================

#worldMaterial = "vacuum"
#rpcMaterial = "vacuum"
#cargoMaterial = "copper"
#length = 100
#width = 100
#thickness = 10
#generation = "monoEnergy_tilt"
#energy = "1GeV"
#angle = "0deg"
#disp = "0cm"
#eventNum = 100000
#
##filename = f"output_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{energy}_{angle}_{disp}_{eventNum}.root"
#filename = f"output_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{eventNum}.root"
##filename = f"output_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{energy}_{angle}_{eventNum}.root"
##filename = "test_option4_limit.root"
#
#utils.checkFileExist(filename)
#
#hits_df, truck_df = utils.returndf(filename)

#nentries = utils.getTotalEventNumber(hits_df)





'''
print("Opening : ", filename )
hits_tree = uproot.open(filename)["hits"]
truck_tree = uproot.open(filename)["truck"]
print("Done opening the file.")
hits_branch = "/eventNumber|muonMomentum|hitTime|hitPosition.|hitPixel.|hitMomentum./"
truck_branch = "/eventNumber|muonMomentum|scatTime|scatPosition.|scatMomentum./"

hits = hits_tree.arrays(filter_name=hits_branch)
truck = truck_tree.arrays(filter_name=truck_branch)

hitPosition_above_list, hitPosition_below_list, truth_scat_list, above_centroid_list, above_dirVec_list, below_centroid_list, below_dirVec_list, poca_list, angle_list, scat_dist_list = utils.returnEventInfo(hits, truck)
'''





# ======================================= Number of scatterings
'''
#outputPDF = f"scattering_number_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{energy}_{angle}_{eventNum}.pdf"  
outputPDF = "test.pdf"
utils.replaceFile(outputPDF)

truth_scat_list = utils.getBranchForAllEvents(truck_df, "scatPositionX_truth")
n_scatterings = [len(truth_scat)-2 for truth_scat in truth_scat_list]

print("Drawing ...")

plt.figure(figsize=(8, 6))
plt.hist(n_scatterings, bins=max(n_scatterings)+2, range=(-0.5, max(n_scatterings)+1.5), edgecolor='black', alpha=0.7)
plt.xlabel("Number of scatterings")
plt.ylabel("Events")
plt.text(0.95, 0.95, f"Entries: {len(n_scatterings)}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
plt.savefig(outputPDF)
'''
'''
for i in range(300):
    filename = f"output_{i}_t0.root"
    hits_df, truck_df = utils.returndf(filename)
    truth_scat_list = utils.getBranchForAllEvents(truck_df, "scatPositionX_truth")
    n_scatterings = [len(truth_scat)-2 for truth_scat in truth_scat_list]
    
    count = n_scatterings.count(0)
    if count > 0:
        print(f"0 scattering exists for {i+1} mm step.")
'''

# ======================================= Visualiztion
'''
output = f"scattering_positions_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{energy}_{angle}_{eventNum}"   
output = "test1"
quantity = 1

pocaResults = utils.poca_reconstruction(hits_df, truck_df)

count = 0
for i in range(len(pocaResults["pocaList"])):

    # Selection
    #if (pocaResults["numScat_list"][i] < 100):
    #    continue    
    ###########

    output_i = output+"_"+str(count)+".html"
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter3d(
        x=[j/10 for j in pocaResults["scatPositionX_truth_list"][i]], 
        y=[j/10 for j in pocaResults["scatPositionZ_truth_list"][i]], 
        z=[j/10 for j in pocaResults["scatPositionY_truth_list"][i]], 
        mode='markers',
        marker=dict(size=1, color='blue', opacity=0.6),
        name="Scattering position"
    ))
    
    fig.add_trace(go.Scatter3d(
        x=[j/10 for j in pocaResults["hitXs_above"][i]], 
        y=[j/10 for j in pocaResults["hitZs_above"][i]], 
        z=[j/10 for j in pocaResults["hitYs_above"][i]], 
        mode='markers',
        marker=dict(size=3, color='red', opacity=0.6),
        name="RPC hits above"
    ))
    
    fig.add_trace(go.Scatter3d(
        x=[j/10 for j in pocaResults["hitXs_below"][i]], 
        y=[j/10 for j in pocaResults["hitZs_below"][i]], 
        z=[j/10 for j in pocaResults["hitYs_below"][i]], 
        mode='markers',
        marker=dict(size=3, color='green', opacity=0.6),
        name="RPC hits below"
    ))
    
    fig.add_trace(go.Scatter3d(
        x=[pocaResults["pocaXList"][i]/10], y=[pocaResults["pocaZList"][i]/10], z=[pocaResults["pocaYList"][i]/10],
        mode='markers',
        marker=dict(size=3, color='black', opacity=1.0),
        name="poca"
    ))
    
    fig.add_trace(utils.create_box((0, 0, 0), (length, width, thickness), color='brown', opacity=0.2))
    #fig.add_trace(create_box((100, 0, 0), (length, 100, thickness), color='brown', opacity=0.2))
    #fig.add_trace(create_box((-100, -150, 0), (length, 100, thickness), color='brown', opacity=0.2))
    #fig.add_trace(create_box((-150, 100, 0), (length, 100, thickness), color='brown', opacity=0.2))

    fig.update_layout(
        title="Interactive 3D Track with Truck and RPC Pixels (Cube View)",
        scene=dict(
            xaxis=dict(title="X (cm)", range=[-length-250, length+250]),
            yaxis=dict(title="Y (cm)", range=[-width-250, width+250]),
            zaxis=dict(title="Z (cm)", range=[-thickness-250, thickness+250]),
            #xaxis=dict(title="X (cm)", range=[-155, 155]),
            #yaxis=dict(title="Y (cm)", range=[-155, 155]),
            #zaxis=dict(title="Z (cm)", range=[-5.5, 5.5]),
            aspectmode="cube"
        ),
        margin=dict(l=0, r=0, t=40, b=0),
    )
    
    above_centroid_swap = [pocaResults["above_centroid_list"][i][0]/10, pocaResults["above_centroid_list"][i][2]/10, pocaResults["above_centroid_list"][i][1]/10]
    below_centroid_swap = [pocaResults["above_centroid_list"][i][0]/10, pocaResults["above_centroid_list"][i][2]/10, pocaResults["above_centroid_list"][i][1]/10]
    above_dirVec_swap = [pocaResults["above_dirVec_list"][i][0], pocaResults["above_dirVec_list"][i][2], pocaResults["above_dirVec_list"][i][1]]
    
    below_dirVec_swap = [pocaResults["below_dirVec_list"][i][0], pocaResults["below_dirVec_list"][i][2], pocaResults["below_dirVec_list"][i][1]]
    
    utils.add_line_auto_range(fig, point=above_centroid_swap, direction=above_dirVec_swap, color='red')
    
    utils.add_line_auto_range(fig, point=below_centroid_swap, direction=below_dirVec_swap, color='green')
    
        
    fig.write_html(output_i)
    
    count += 1
    
    if count >= quantity:
        break
'''

# ======================================= 2D Visualization
#
#thickness = 1
#
#filename = "test0.root"
#utils.checkFileExist(filename)
#hits_df, truck_df = utils.returndf(filename)
#nentries = utils.getTotalEventNumber(hits_df)
#
#filename_limit = "test1.root"
#utils.checkFileExist(filename_limit)
#hits_df_limit, truck_df_limit = utils.returndf(filename_limit)
#
#
#output = "2D_scattering_positions_test0_test1"
#
#quantity = 5
#
#pocaResults = utils.poca_reconstruction(hits_df, truck_df)
#pocaResults_limit = utils.poca_reconstruction(hits_df_limit, truck_df_limit)
#
#count = 0


#for i in range(len(pocaResults["pocaList"])):
#
#    output_i = output+"_"+str(count)+".pdf"
#    
#    x = [j/10 for j in pocaResults["scatPositionX_truth_list"][i]]
#    z = [j/10 for j in pocaResults["scatPositionY_truth_list"][i]]
#    plt.figure(figsize=(8, 5))
#    plt.scatter(x, z, c='blue', label='Scattering point')
#    plt.plot(x, z, 'r--', alpha=0.5) 
#    plt.xlim(-0.002, 0.002)
#    plt.xlabel('x-axis', fontsize=12)
#    plt.ylabel('z-axis', fontsize=12)
#    plt.grid(True, linestyle='--', alpha=0.7)
#    
#    plt.savefig(output_i)
#    plt.close()
#    
#    count += 1
#    
#    if count >= quantity:
#        break


#for i in range(len(pocaResults["pocaList"])):
#
#    output_i = output+"_"+str(count)+".pdf"
#
#    
#    
#    x = [j/10 for j in pocaResults["scatPositionX_truth_list"][i]]
#    x.append(pocaResults["hitXs_below"][i][0] / 10)
#    z = [j/10 for j in pocaResults["scatPositionY_truth_list"][i]]
#    z.append(pocaResults["hitYs_below"][i][0] / 10)
#    
#    x_limit = [j/10 for j in pocaResults_limit["scatPositionX_truth_list"][i]]
#    x_limit.append(pocaResults_limit["hitXs_below"][i][0] / 10)
#    z_limit = [j/10 for j in pocaResults_limit["scatPositionY_truth_list"][i]]
#    z_limit.append(pocaResults_limit["hitYs_below"][i][0] / 10)
#    
#    plt.figure(figsize=(8, 5))
#    plt.scatter(x, z, c='blue', marker='o', alpha=0.7, label='No step limit')
#    plt.scatter(x_limit, z_limit, c='red', marker='x', alpha=0.7, label='With step limit')
#    plt.plot(x, z, '--', color='royalblue', alpha=0.5)  
#    plt.plot(x_limit, z_limit, ':', color='crimson', alpha=0.5)  
#    plt.xlim(-0.002, 0.002)
#    plt.ylim(-0.51, 0.51)
#    plt.xlabel('x-axis', fontsize=12)
#    plt.ylabel('z-axis', fontsize=12)
#    plt.grid(True, linestyle='--', alpha=0.7)
#    plt.legend()
#    plt.savefig(output_i)
#    plt.close()
#    
#    count += 1
#    
#    if count >= quantity:
#        break


# ======================================= POCA point distribution
'''
for cnt in range(5):

    filename = f"output_{cnt}_copper.root"
    
    utils.checkFileExist(filename)
    
    hits_df, truck_df = utils.returndf(filename)
    
    nentries = utils.getTotalEventNumber(hits_df)
    
    length = 100
    width = 100
    thickness = 100
    
    #outputPDF = f"poca_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{eventNum}.pdf" 
    outputPDF = f"poca_{cnt}_copper.pdf"
    utils.replaceFile(outputPDF)
    
    pocaResults = utils.poca_reconstruction(hits_df, truck_df)
    
    pocaX_swapped_cm = [i/10 for i in pocaResults["pocaXList"]]
    pocaY_swapped_cm = [i/10 for i in pocaResults["pocaZList"]]
    pocaZ_swapped_cm = [i/10 for i in pocaResults["pocaYList"]]
    
    space = 5 
    
    with PdfPages(outputPDF) as pdf:    
    
        plt.figure(figsize=(8,6))
        plt.hist(pocaX_swapped_cm, bins=200, range=(-length-space, length+space), edgecolor='black', alpha=0.7)
        plt.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
        plt.xlabel("Scattering x [cm]")
        plt.ylabel("Events")
        plt.text(0.95, 0.95, f"Entries: {len(pocaX_swapped_cm)}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        plt.text(0.95, 0.88, f"Mean: {np.mean(pocaX_swapped_cm):.4f}\nStd Dev: {np.std(pocaX_swapped_cm):.3f}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        pdf.savefig()
        plt.close()
    
        plt.figure(figsize=(8,6))
        plt.hist(pocaY_swapped_cm, bins=200, range=(-width-space, width+space), edgecolor='black', alpha=0.7)
        plt.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
        plt.xlabel("Scattering y [cm]")
        plt.ylabel("Events")
        plt.text(0.95, 0.95, f"Entries: {len(pocaY_swapped_cm)}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        plt.text(0.95, 0.88, f"Mean: {np.mean(pocaY_swapped_cm):.4f}\nStd Dev: {np.std(pocaY_swapped_cm):.3f}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        pdf.savefig()
        plt.close()
    
        plt.figure(figsize=(8,6))
        plt.hist(pocaZ_swapped_cm, bins=200, range=(-thickness-space, thickness+space), edgecolor='black', alpha=0.7)
        plt.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
        plt.xlabel("Scattering z [cm]")
        plt.ylabel("Events")
        plt.text(0.95, 0.95, f"Entries: {len(pocaZ_swapped_cm)}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        plt.text(0.95, 0.88, f"Mean: {np.mean(pocaZ_swapped_cm):.4f}\nStd Dev: {np.std(pocaZ_swapped_cm):.3f}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        pdf.savefig()
        plt.close()
'''

'''
def gaussian(x, mean, std, amplitude):
    return amplitude * np.exp(-(x - mean)**2 / (2 * std**2))

def crystal_ball(x, alpha_l, n_l, alpha_r, n_r, mean, sigma, N):

    x = (x - mean) / sigma 

    A_l = (n_l / np.abs(alpha_l))**n_l * np.exp(-0.5 * alpha_l**2)
    B_l = n_l / np.abs(alpha_l) - np.abs(alpha_l)
    left_tail = A_l * (B_l - x)**(-n_l)

    A_r = (n_r / np.abs(alpha_r))**n_r * np.exp(-0.5 * alpha_r**2)
    B_r = n_r / np.abs(alpha_r) - np.abs(alpha_r)
    right_tail = A_r * (B_r + x)**(-n_r)

    gaussian_core = np.exp(-0.5 * x**2)

    result = np.where(x < -alpha_l, left_tail,
                     np.where(x > alpha_r, right_tail, gaussian_core))

    return N * result

matList = ["copper", "water", "tungsten"]

for mat in matList:

    energies = ['1GeV', '5GeV', '10GeV']
    colors = ['#1f77b4', '#2ca02c', '#d62728']
    
    all_parameters = {}
    for energy in energies:
        all_parameters[energy] = {
            'steps': [],
            'mean': [], 'mean_err': [],
            'sigma': [], 'sigma_err': [],
            'alpha_l': [], 'alpha_l_err': [],
            'n_l': [], 'n_l_err': [],
            'alpha_r': [], 'alpha_r_err': [],
            'n_r': [], 'n_r_err': [],
            'amp': [], 'amp_err': []
        }
    
    outputPDF = f"{mat}_1_5_10_GeV_poca_Z_fit_DCB.pdf"
    utils.replaceFile(outputPDF)
    
    stepList = []
    for i in range(1, 20, 1):
        stepList.append(i / 10) 
    
    
    with PdfPages(outputPDF) as pdf:    
        
        for energy in energies:
        
            for i in range(0,19):
                
                step_mm = stepList[i]
            
                print(f"Processing file {i} ...")
                
                filename = f"stepLimitValidation/output_{i}_{mat}_{energy}.root"
                        
                hits_df, truck_df = utils.returndf(filename)
                
                pocaResults = utils.poca_reconstruction(hits_df, truck_df)
                pocaZ_swapped_cm = [i/10 for i in pocaResults["pocaYList"]] 
            
                plt.figure(figsize=(8,6))
                
                y_data, x_edges, _ = plt.hist(pocaZ_swapped_cm, bins=200, range=(-1.5, 1.5), edgecolor='black', alpha=0.7, label="Data")
                bin_centers = (x_edges[:-1] + x_edges[1:]) / 2
                bin_width = x_edges[1] - x_edges[0]
                #initial_guess = [np.mean(pocaZ_swapped_cm), np.std(pocaZ_swapped_cm), max(y_data)]
                initial_guess = [0.5, 4, 0.5, 4, np.mean(pocaZ_swapped_cm), np.std(pocaZ_swapped_cm), max(y_data)]
                try:
                    #popt, pcov = curve_fit(gaussian, bin_centers, y_data, p0=initial_guess)
                    #mean_fit, std_fit, amp_fit = popt
                    
                    lower_bounds = [-np.inf, 1, -np.inf, 1, -np.inf, 0, 0] 
                    upper_bounds = [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]
                    popt, pcov = curve_fit(crystal_ball, bin_centers, y_data, p0=initial_guess, bounds=(lower_bounds, upper_bounds))
                    perr = np.sqrt(np.diag(pcov))
                    alpha_l, n_l, alpha_r, n_r, mean_fit, sigma_fit, amp_fit = popt
                    alpha_l_err, n_l_err, alpha_r_err, n_r_err, mean_err, sigma_err, amp_err = perr
                    x_fit = np.linspace(-1.5, 1.5, 1000)
                    #y_fit = gaussian(x_fit, *popt)
                    y_fit = crystal_ball(x_fit, *popt)
                    #plt.plot(x_fit, y_fit, 'r-', lw=2, label=f'Gaussian Fit\nMean: {mean_fit:.4f}\nStd Dev: {std_fit:.3f}')
                    fit_label = (
                        f'DCB Fit\n'
                        f'Mean: {mean_fit:.4f} +- {mean_err:.4f}\n'
                        f'Sigma: {sigma_fit:.3f} +- {sigma_err:.3f}\n'
                        f'alpha_l: {alpha_l:.3f} +- {alpha_l_err:.3f}\n'
                        f'n_l: {n_l:.3f} +- {n_l_err:.3f}\n' 
                        f'alpha_r: {alpha_r:.3f} +- {alpha_r_err:.3f}\n'
                        f'n_r: {n_r:.3f} +- {n_r_err:.3f}'
                    )
                    
                    plt.plot(x_fit, y_fit, 'r-', lw=2, label=fit_label)
                    plt.legend(loc='upper left', fontsize=10)
                    
                    
                    all_parameters[energy]['steps'].append(step_mm)
                    all_parameters[energy]['mean'].append(mean_fit)
                    all_parameters[energy]['mean_err'].append(mean_err)
                    all_parameters[energy]['sigma'].append(sigma_fit)
                    all_parameters[energy]['sigma_err'].append(sigma_err)
                    all_parameters[energy]['alpha_l'].append(alpha_l)
                    all_parameters[energy]['alpha_l_err'].append(alpha_l_err)
                    all_parameters[energy]['n_l'].append(n_l)
                    all_parameters[energy]['n_l_err'].append(n_l_err)
                    all_parameters[energy]['alpha_r'].append(alpha_r)
                    all_parameters[energy]['alpha_r_err'].append(alpha_r_err)
                    all_parameters[energy]['n_r'].append(n_r)
                    all_parameters[energy]['n_r_err'].append(n_r_err)
                    all_parameters[energy]['amp'].append(amp_fit)
                    all_parameters[energy]['amp_err'].append(amp_err)
                    
                    
                except RuntimeError:
                    plt.text(0.05, 0.95, "Fit failed", transform=plt.gca().transAxes, fontsize=12, color='red')
                    all_parameters[energy]['steps'].append(step_mm)
                    for key in ['mean', 'sigma', 'alpha_l', 'n_l', 'alpha_r', 'n_r', 'amp']:
                        all_parameters[energy][key].append(np.nan)
                        all_parameters[energy][f'{key}_err'].append(np.nan)
                
                
                #plt.hist(pocaZ_swapped_cm, bins=200, range=(-1.5, 1.5), edgecolor='black', alpha=0.7)
                plt.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
                plt.xlabel("Scattering z [cm]")
                plt.ylabel("Events")
                plt.text(0.95, 0.95, f"Step limit: {round(step_mm, 2)} mm", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
                plt.text(0.95, 0.88, f"Entries: {len(pocaZ_swapped_cm)}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
                #plt.text(0.95, 0.81, f"Mean: {np.mean(pocaZ_swapped_cm):.4f}\nStd Dev: {np.std(pocaZ_swapped_cm):.3f}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
                pdf.savefig()
                plt.close()
            
    
    def plot_parameter_trend(param_name, ylabel, ax):
        for energy, color in zip(energies, colors):
            params = all_parameters[energy]
            ax.errorbar(
                x=params['steps'],
                y=params[param_name],
                yerr=params[f'{param_name}_err'],
                fmt='o', color=color,
                markersize=6, capsize=4, capthick=2,
                label=f"{energy}"                
            )
            ax.legend(loc='best')
            ax.set_xlabel('Step Limit (mm)', fontsize=12)
            ax.set_ylabel(ylabel, fontsize=12)
            ax.grid(True, alpha=0.3)
            if param_name == "mean":
                ax.set_ylim(-0.5, 0.5)
            ax.set_xlim(0, 2.1)
            ax.axhline(0, color='gray', linestyle='--', linewidth=1)
    
    fig = plt.figure(figsize=(16, 20))
    gs = gridspec.GridSpec(4, 2, hspace=0.4, wspace=0.3)
    
    params_config = [
        ('mean', 'Mean Position [cm]'),
        ('sigma', 'Sigma [cm]'),
        ('alpha_l', 'Left Alpha'),
        ('n_l', 'Left n'),
        ('alpha_r', 'Right Alpha'), 
        ('n_r', 'Right n'),
        ('amp', 'Amplitude')
    ]
    
    for idx, (param, label) in enumerate(params_config):
        ax = fig.add_subplot(gs[idx])
        plot_parameter_trend(param, label, ax)
    
    summary_pdf = f"{mat}_fit_parameters_summary.pdf"
    plt.savefig(summary_pdf, bbox_inches='tight')
    plt.close()

'''


# ======================================= Scattering angle distribution
'''
outputPDF = f"scatAngle_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{eventNum}.pdf" 


utils.replaceFile(outputPDF)

#pocaResults = utils.poca_reconstruction(hits_df, truck_df)

proj_angle_list = pocaResults["angle_list"]
#proj_angle_list = [i for i, j in zip(pocaResults["proj_angle_list"], pocaResults["numScat_list"]) if j == 0]

with PdfPages(outputPDF) as pdf:
    plt.figure(figsize=(8,6))
    plt.hist(proj_angle_list, bins=200, range=(-0.5, 50), edgecolor='black', alpha=0.7)
    plt.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
    plt.xlabel("Scattering angle [mrad]")
    plt.ylabel("Events")
    plt.text(0.95, 0.95, f"Entries: {len(proj_angle_list)}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
    plt.text(0.95, 0.88, f"Mean: {np.mean(proj_angle_list):.4f}\nStd Dev: {np.std(proj_angle_list):.3f}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
    pdf.savefig()
    plt.close()
'''
# ====================================== Scattering distance distribution

for cnt in range(5):

    filename = f"output_{cnt}_copper.root"
    
    utils.checkFileExist(filename)
    
    hits_df, truck_df = utils.returndf(filename)
    
    nentries = utils.getTotalEventNumber(hits_df)
    
    length = 100
    width = 100
    thickness = 100
    

    #outputPDF = f"scatDist_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{eventNum}.pdf" 
    outputPDF = f"scatDist_{cnt}_copper_0p1.pdf"  
    
    utils.replaceFile(outputPDF)
    
    pocaResults = utils.poca_reconstruction(hits_df, truck_df)
    
    scat_dist_list = pocaResults["scat_dist"]
    
    with PdfPages(outputPDF) as pdf:
        plt.figure(figsize=(8,6))
        plt.hist(scat_dist_list, bins=200, range=(0, 0.01), edgecolor='black', alpha=0.7, log=True)
        plt.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
        plt.xlabel("Scattering distance [cm]")
        plt.ylabel("Events")
        plt.text(0.95, 0.95, f"Entries: {len(scat_dist_list)}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        plt.text(0.95, 0.88, f"Mean: {np.mean(scat_dist_list):.4f}\nStd Dev: {np.std(scat_dist_list):.3f}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        pdf.savefig()
        plt.close()



# ======================================= POCA point Z distribution with fixed number of scattering
'''
outputPDF = f"pocaZ_separate_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{energy}_{angle}_{disp}_{eventNum}.pdf"  
#outputPDF = "test.pdf"

utils.replaceFile(outputPDF)

with PdfPages(outputPDF) as pdf:
    
    for numScat in range(50):
        
        filtered_poca_Z = [sublist[-1] for i, sublist in enumerate(poca_list) if (len(truth_scat_list[i])-2 == numScat)]
        
        if len(filtered_poca_Z) <= 1:
            continue
            
        plt.figure(figsize=(8,6))
        plt.hist(filtered_poca_Z, bins=200, range=(-thickness-15, thickness+15), edgecolor='black', alpha=0.7)
        plt.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
        plt.xlabel("Scattering z [cm]")
        plt.ylabel("Events")
        plt.text(0.95, 0.95, f"Entries: {len(filtered_poca_Z)}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        plt.text(0.95, 0.88, f"Mean: {np.mean(filtered_poca_Z):.4f}\nStd Dev: {np.std(filtered_poca_Z):.3f}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        plt.text(0.95, 0.80, f"Number of scattering: {numScat}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
        pdf.savefig()
        plt.close()
'''

# ======================================= POCA point Z distribution against displacement
'''
worldMaterial = "vacuum"
rpcMaterial = "vacuum"
cargoMaterial = "copper"
length = 100
width = 100
thickness = 100
generation = "monoEnergy_tilt"
energy = "10GeV"
angle = "30deg"
disps = ["10cm", "30cm", "50cm", "70cm", "90cm", "110cm", "130cm", "150cm"]
eventNum = 50000

outputPDF = f"pocaZ_inclined_filter0Scattering_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{energy}_{angle}_{eventNum}.pdf"  
#outputPDF = "test.pdf"

utils.replaceFile(outputPDF)

with PdfPages(outputPDF) as pdf:

    plt.figure(figsize=(10, 6))
    colors = ['red', 'orange', 'yellow', 'green', 'blue', 'purple', 'black', 'magenta']
    labels = ["10cm", "30cm", "50cm", "70cm", "90cm", "110cm", "130cm", "150cm"]
    #boundary = [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [10, 100]]
    
    all_hist_data = []

    for i in range(len(disps)):
        
        filename = f"output_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{energy}_{angle}_{disps[i]}_{eventNum}.root"
        
        file_exists = os.path.exists(filename)
        if not file_exists:
            print(f"File '{filename}' does not exist.")
            sys.exit(1)
        
        print("Opening : ", filename )
        hits_tree = uproot.open(filename)["hits"]
        truck_tree = uproot.open(filename)["truck"]
        print("Done opening the file.")
        hits_branch = "/eventNumber|muonMomentum|hitTime|hitPosition.|hitPixel.|hitMomentum./"
        truck_branch = "/eventNumber|muonMomentum|scatTime|scatPosition.|scatMomentum./"
        
        hits = hits_tree.arrays(filter_name=hits_branch)
        truck = truck_tree.arrays(filter_name=truck_branch)
        
        hitPosition_above_list, hitPosition_below_list, truth_scat_list, above_centroid_list, above_dirVec_list, below_centroid_list, below_dirVec_list, poca_list, angle_list = utils.returnEventInfo(hits, truck)
        
        filtered_poca_Z = [sublist[-1] for i, sublist in enumerate(poca_list) if (len(truth_scat_list[i])-2 > 0)]                    

        counts, bins, _ = plt.hist(
            filtered_poca_Z, bins=200, range=(-thickness-10, thickness+10),
            color=colors[i], alpha=0.7, label=labels[i],
            edgecolor='black', linewidth=0.5,
            stacked=True
        )
        all_hist_data.append(counts)
        
        
    plt.xlabel("Scattering z [cm]")
    plt.ylabel("Events")
    plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.3)
    
    #if all_hist_data:
    #    max_count = max(np.sum(counts) for counts in zip(*all_hist_data))
    #    plt.ylim(0, max_count * 1.1)
    plt.ylim(0, 10000)
    
    plt.tight_layout()
    pdf.savefig()
    plt.close()
'''



# ======================================= POCA point Z distribution with fixed number of scattering
'''
outputPDF = f"pocaZ_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{energy}_{angle}_{disp}_{eventNum}.pdf"  
#outputPDF = "test.pdf"

utils.replaceFile(outputPDF)

with PdfPages(outputPDF) as pdf:

    plt.figure(figsize=(10, 6))
    colors = ['red', 'orange', 'yellow', 'green', 'blue', 'purple']
    labels = ['1-2 Scatterings', '3-4 Scattering', '5-6 Scatterings', '7-8 Scattering', '9-10 Scatterings', '>= 11 Scatterings']
    boundary = [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [10, 100]]
    
    all_hist_data = []

    for numScat in range(len(boundary)):

        filtered_poca_Z = [sublist[-1] for i, sublist in enumerate(poca_list) if (len(truth_scat_list[i])-2 >= boundary[numScat][0] and len(truth_scat_list[i])-2 <= boundary[numScat][1])]
        
        if len(filtered_poca_Z) <= 1:
            continue
        
        counts, bins, _ = plt.hist(
            filtered_poca_Z, bins=200, range=(-thickness-10, thickness+10),
            color=colors[numScat], alpha=0.7, #label=labels[numScat],
            edgecolor='black', linewidth=0.5,
            stacked=True
        )
        all_hist_data.append(counts)
        
        
    plt.xlabel("Scattering z [cm]")
    plt.ylabel("Events")
    #plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.3)
    
    #if all_hist_data:
    #    max_count = max(np.sum(counts) for counts in zip(*all_hist_data))
    #    plt.ylim(0, max_count * 1.1)
    plt.ylim(0, 15000)
        
    for numScat in range(len(boundary)):
        filtered_poca_Z = [
            sublist[-1] for i, sublist in enumerate(poca_list) 
            if (len(truth_scat_list[i])-2 >= boundary[numScat][0] and len(truth_scat_list[i])-2 <= boundary[numScat][1])
        ]
        if len(filtered_poca_Z) > 1:
            plt.text(
                0.95, 0.85 - numScat*0.1, 
                f"{labels[numScat]}: \nMean = {np.mean(filtered_poca_Z):.4f} cm\nStd Dev = {np.std(filtered_poca_Z):.4f} cm",
                ha='right', va='top', transform=plt.gca().transAxes, 
                fontsize=10, color=colors[numScat]
            )
    
    plt.tight_layout()
    pdf.savefig()
    plt.close()
'''
# ======================================= POCA z vs Number of scatterings
'''
outputPDF = f"pocaZ_vs_numScat_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{energy}_{angle}_{eventNum}.pdf"  
#outputPDF = f"test.pdf"

utils.replaceFile(outputPDF)

numScat = [len(sublist) - 2 for sublist in truth_scat_list]
pocaZ = [sublist[2] for sublist in poca_list]

groups = {}
for i, val in zip(numScat, pocaZ):
    if i not in groups:
        groups[i] = []
    groups[i].append(val)

x_values = sorted(groups.keys())
means = [np.mean(groups[x]) for x in x_values]

stds = []
std_errors = []
for x in x_values:
    samples = groups[x]
    n = len(samples)
    if n > 1:
        std = np.std(samples, ddof=1)
        std_err = std / np.sqrt(n)
    else:
        std = 0
        std_err = 0
    stds.append(std)
    std_errors.append(std_err)

plt.figure(figsize=(8, 5))
plt.plot(x_values, means, 'o:', markersize=8)
mask = np.array([len(groups[x]) > 1 for x in x_values])
if any(mask):
    plt.errorbar(
        np.array(x_values)[mask],
        np.array(means)[mask],
        yerr=np.array(std_errors)[mask],
        fmt='none',
        capsize=5,
        ecolor='black',
    )

plt.xlabel('Number of scattering', fontsize=12)
plt.ylabel('Average POCA z [cm]', fontsize=12)
plt.xticks(np.arange(0, max(x_values)+1))
plt.grid(True, linestyle='--', alpha=0.5)
plt.savefig(outputPDF)
'''

# ======================================= Scattering distance vs Number of scatterings
'''
outputPDF = f"scatDist_vs_numScat_{worldMaterial}_{rpcMaterial}_{cargoMaterial}_{length}cm_{width}cm_{thickness}cm_{generation}_{energy}_{angle}_{eventNum}.pdf"  
#outputPDF = f"test.pdf"

utils.replaceFile(outputPDF)

numScat = [len(sublist) - 2 for sublist in truth_scat_list]

groups = {}
for i, val in zip(numScat, scat_dist_list):
    if i not in groups:
        groups[i] = []
    groups[i].append(val)

x_values = sorted(groups.keys())
means = [np.mean(groups[x]) for x in x_values]

stds = []
std_errors = []
for x in x_values:
    samples = groups[x]
    n = len(samples)
    if n > 1:
        std = np.std(samples, ddof=1)
        std_err = std / np.sqrt(n)
    else:
        std = 0
        std_err = 0
    stds.append(std)
    std_errors.append(std_err)

plt.figure(figsize=(8, 5))
plt.plot(x_values, means, 'o:', markersize=8)
mask = np.array([len(groups[x]) > 1 for x in x_values])
if any(mask):
    plt.errorbar(
        np.array(x_values)[mask],
        np.array(means)[mask],
        yerr=np.array(std_errors)[mask],
        fmt='none',
        capsize=5,
        ecolor='black',
    )

plt.xlabel('Number of scattering', fontsize=12)
plt.ylabel('Average scattering distance [cm]', fontsize=12)
plt.xticks(np.arange(0, max(x_values)+1))
plt.grid(True, linestyle='--', alpha=0.5)
plt.savefig(outputPDF)
'''
