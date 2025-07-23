# author: Lisa Bast
# date: 2025-05-20,  09:43:12
# version: 0.0.1
# about: clustermap median AUC value per donor and significant regulon, 
#        clustermap median AUC value per donor for each cell type
#AUC per regulon per cell = enrichment score for the regulon's target genes = approximation for the activity of the regulon

import sys
import os
import pandas as pd
import seaborn as sns
import numpy as np

path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

# paths and files:
path_main = path_code.replace("5e_regulatory_network_inference\\script","")
path_main_project = path_code.replace("script","")
path_AUC_results = path_main_project + "output/pyscenic/"
path_figures = path_AUC_results+"figures/AUC_plots/"

results_file_sign = "median_regulon_AUC_donor_celltype_sign.csv"

##clustermap median AUC value per donor and significant regulon
#read data frames
DF_AUC_pivot,DF_cond_donor = ut.get_sign_AUC_pivot_and_cond_donor(path_AUC_results,results_file_sign)

DF_AUC_pivot_no_nan = DF_AUC_pivot.fillna(0)
DF_AUC_pivot_no_nan = pd.merge(left=DF_AUC_pivot_no_nan,right=DF_cond_donor,left_index=True,right_index=True,how="left")

#define colors
row_colors,groups = ut.get_cond_colors(DF_AUC_pivot_no_nan)
cmap_light_tur = sns.color_palette("light:#5A9", as_cmap=True)

ut.plot_and_save_clustermap(DF_AUC_pivot_no_nan,row_colors,cmap_light_tur, 'median AUC', 'Significant Regulons', "median_AUC_sign_all_celltypes",path_figures)

##z-score per regulon across donors
g = ut.plot_and_save_clustermap(DF_AUC_pivot_no_nan,row_colors,"vlag_r", 'Z-score of median AUC per donor', 'Significant Regulons', "z_score_median_AUC_sign_all_celltypes",path_figures)

##resort the y axis into groups and plot again
samples_ordered = ut.resort_yaxis_labels_by_group(g,groups)
DF_AUC_pivot_no_nan_reordered = DF_AUC_pivot_no_nan.reindex(samples_ordered)
_,DF_cond_donor = ut.get_sign_AUC_pivot_and_cond_donor(path_AUC_results,results_file_sign)
DF_AUC_pivot_no_nan_reordered = pd.merge(left=DF_AUC_pivot_no_nan_reordered,right=DF_cond_donor,left_index=True,right_index=True,how="left")

row_colors,_ = ut.get_cond_colors(DF_AUC_pivot_no_nan_reordered)

ut.plot_and_save_clustermap(DF_AUC_pivot_no_nan_reordered,row_colors,"vlag_r", 'z-score_per_regulon_of_sign_median_AUC_all_celltypes_within_group_clustering', 'Significant Regulons', 'z-score_per_regulon_of_sign_median_AUC_all_celltypes_within_group_clustering',path_figures)


##clustermap median AUC value per donor and regulon for each cell type
DF_AUC_all = ut.get_AUC_pivot(path_AUC_results,"median_regulon_AUC_donor_celltype.csv")
_,DF_cond_donor = ut.get_sign_AUC_pivot_and_cond_donor(path_AUC_results,results_file_sign)
CTs = np.unique(DF_AUC_all["celltype"].tolist())
for ct in CTs:
    #only keep current ct
    DF_AUC_ct = DF_AUC_all[DF_AUC_all["celltype"]==ct].copy()
    DF_AUC_pivot = pd.pivot_table(DF_AUC_ct[["regulon","median_AUC","Donor","condition"]],columns="regulon",values="median_AUC",index=["Donor"])
    DF_AUC_pivot_no_nan = DF_AUC_pivot.fillna(0)
    #remove regulons that are zero for all donors:
    DF_AUC_pivot_sparse = DF_AUC_pivot_no_nan[DF_AUC_pivot_no_nan.columns[DF_AUC_pivot_no_nan.sum()!=0]]
    DF_AUC_pivot_sparse = pd.merge(left=DF_AUC_pivot_sparse,right=DF_cond_donor,left_index=True,right_index=True,how="left")

    row_colors,groups = ut.get_cond_colors(DF_AUC_pivot_sparse)

    cmap_light_tur = sns.color_palette("light:#5A9", as_cmap=True)

    ut.plot_and_save_clustermap(DF_AUC_pivot_sparse,row_colors,cmap_light_tur,'median_AUC_'+ct, 'All Regulons '+ct, 'median_AUC_'+ct, path_figures)
