# author: Lisa Bast
# date: 2025-05-20,  16:35:07
# version: 0.0.1
# about: utils

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_sign_AUC_pivot_and_cond_donor(path_AUC_results,results_file_sign):
    DF_AUC = pd.read_csv(path_AUC_results+results_file_sign)
    DF_AUC["Celltype_regulon"] = DF_AUC["celltype"]+' '+DF_AUC["regulon"]
    DF_AUC_pivot = pd.pivot_table(DF_AUC[["Celltype_regulon","median_AUC","Donor","condition"]],columns="Celltype_regulon",values="median_AUC",index=["Donor"])

    DF_cond_donor= DF_AUC[["Donor","condition"]]
    DF_cond_donor.drop_duplicates(inplace=True)
    DF_cond_donor.set_index("Donor",inplace=True)

    return DF_AUC_pivot,DF_cond_donor

def get_AUC_pivot(path_AUC_results,results_file):
    DF_AUC_all = pd.read_csv(path_AUC_results+results_file)

    return DF_AUC_all

def get_cond_colors(DF):
    groups = DF.pop("condition")
    lut = {"CTRL":"steelblue","SCZ":"darkorange"}#dict(zip(groups.unique(), "rbg"))
    row_colors = groups.map(lut)
    return row_colors,groups


#plot clustermap
def plot_and_save_clustermap(DF,row_colors,cmap,cbar_label, title, option, path_figures):
    if option == "median_AUC_sign_all_celltypes":
        g = sns.clustermap(DF,row_colors=row_colors, cmap = cmap,figsize=(7,17),cbar_kws={'label': cbar_label}) #z_score=0,cmap="vlag",center=0,row_colors=row_colors)
    elif option=='z-score_per_regulon_of_sign_median_AUC_all_celltypes':
        g = sns.clustermap(DF,row_colors=row_colors, zscore=1, center=0, cmap = cmap,figsize=(7,17),cbar_kws={'label': cbar_label}) #z_score=0,cmap="vlag",center=0,row_colors=row_colors)
    elif option=='z-score_per_regulon_of_sign_median_AUC_all_celltypes_within_group_clustering':
        g = sns.clustermap(DF, row_colors=row_colors, z_score=1, center=0, row_cluster=False, cmap = cmap,figsize=(7,17),cbar_kws={'label': cbar_label}) #z_score=0,cmap="vlag",center=0,row_colors=row_colors)
    else:
        g = sns.clustermap(DF,row_colors=row_colors, cmap = cmap,figsize=(10,17),cbar_kws={'label': cbar_label})
    #highlight the significant regulons
    g.figure.suptitle(title) 
    g.ax_cbar.set_position((0.01, .8, .03, .1))

    plt.savefig(path_figures+option+'.pdf', dpi=400)
    return g


##resort the y axis into groups and plot again
def resort_yaxis_labels_by_group(g,groups):
    ylabels = [t.get_text() for t in g.ax_heatmap.yaxis.get_majorticklabels()]
    ctrl_samples_ordered = []
    scz_samples_ordered = []
    for l in ylabels:
        if groups[l]=="CTRL":
            ctrl_samples_ordered = ctrl_samples_ordered + [l]
        else:
            scz_samples_ordered = scz_samples_ordered + [l]
    samples_ordered = ctrl_samples_ordered + scz_samples_ordered
    return samples_ordered
