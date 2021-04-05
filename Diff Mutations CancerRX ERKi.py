import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os
import math

"""Set work directory"""
os.chdir('C:\\Lab\\Python\\Github\\ERK inhibitors')
"""Data file listing cell types"""
df_cell_types= pd.read_excel('CancerRX_Cell_lines.xlsx')
"""File containing genetic features data (data obtained from CancerRX database)"""
df_genetic= pd.read_excel('PANCANCER_Genetic_feature__Thu Nov  5 22_20_38 2020.xlsx')
mutations= list(set(dict.fromkeys(df_genetic['genetic_feature'])))
"""Data set containing drugs' IC50 values from CancerRX database"""
drugs= {'Ulix': 'Ulix_PANCANCER.xlsx',
        'SCH': 'SCH_PANCANCER.xlsx',
        'VX11e': 'VX11e_PANCANCER.xlsx'}
"""Specify which cell type, drug, and gene set to use in analysis"""
cell_type= 'Solid' # Solid or Blood
drug= 'VX11e' # Ulix, SCH or VX11e

df_cell_types= pd.read_excel('CancerRX_Cell_lines.xlsx')
df_drug= pd.read_excel(drugs[drug])
cells= set(df_drug['Cell line name']) & set(df_cell_types[cell_type]) & set(df_genetic['cell_line_name']) #selects cell line with genetic features and IC50 data available
df_drug['IC50_exp']= df_drug['IC50'].apply(lambda x: math.exp(x)) #converts IC50 values from ln() to linear values to avoid negative values
df_drug.set_index('Cell line name', inplace=True)
p_list=list()
dif_list=list()
"""For each genetic feature (feat) cell lines are divided in mutant (mut) and non-mutant (wt),
    and then IC50 values were compared between these groups. p-value calculated for each feature
    were then adjusted by FDR."""
for feat in mutations:
    df_feat= df_genetic[df_genetic['genetic_feature']==feat]
    df_mut= df_feat[df_feat['is_mutated']==1]
    df_wt= df_feat[df_feat['is_mutated']==0]
    cells_mut= list(set(dict.fromkeys(df_mut['cell_line_name'])) & cells)
    cells_wt= list(set(dict.fromkeys(df_wt['cell_line_name'])) & cells)
    stat, p= mannwhitneyu(df_drug.loc[cells_mut]['IC50_exp'], df_drug.loc[cells_wt]['IC50_exp'], alternative='two-sided')
    dif= df_drug.loc[cells_mut]['IC50_exp'].mean() / df_drug.loc[cells_wt]['IC50_exp'].mean()
    p_list.append(p)
    dif_list.append(dif)
p_adj= multipletests(pvals=p_list, alpha=0.01, method='fdr_tsbky')
df_out= pd.DataFrame({'Feature': mutations, 'p-value': p_list, 'p-adj': p_adj[1],
                      'diffference mut vs wt': dif_list, 'significant': p_adj[0]})

out_dir= os.path.join(os.getcwd(), 'Differential Analysis\\')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
df_out.to_excel(out_dir+'Diff_mutations_FDR_'+drug+'_'+cell_type+'.xlsx')
