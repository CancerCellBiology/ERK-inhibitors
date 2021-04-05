import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV
import os
from sklearn.linear_model import ElasticNet
from scipy.stats import zscore
import math

def elastic_param (X, y, cross_val):
    """Computes elastic net parameters and returns best parameters
    X- arguments matrix; y- estimated value; cross_val- number of folds for cross validation"""
    elastic=ElasticNet(normalize=True, tol=tolerance)
    search=GridSearchCV(estimator=elastic,param_grid={'alpha':[0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 1],'l1_ratio':[.1,.2,.4,.6,.8]},scoring='neg_mean_squared_error',n_jobs=1,refit=True,cv=cross_val)
    search.fit(X,y)
    print(search.best_params_)
    print(abs(search.best_score_))
    return search.best_params_

def perform_elastic (X, y, a, l1, n):
    """Performs elastic net with specified parameters and 
    returns a DataFrame with mean elastic net coefficients
    for all cycles and z-scores for each argument.
    X- arguments matrix; y- estimated value; a- alpha penalty; l1- L1 and L2 combined penalty;
    n- number of cycles elastic net will be performed"""
    elastic=ElasticNet(normalize=True ,alpha=a, l1_ratio=l1, tol=tolerance, selection='random')
    df_elastic=pd.DataFrame()
    for x in range(0, n):
        elastic.fit(X,y)
        coef_dict_baseline = {}
        for coef, feat in zip(elastic.coef_,X.columns):
            coef_dict_baseline[feat] = coef
        df_temp = pd.DataFrame.from_dict(coef_dict_baseline, orient='index')
        df_elastic = pd.concat([df_elastic, df_temp], axis=1)
    df_out = pd.DataFrame(index=df_elastic.index)
    df_out['Result'] = df_elastic.mean(axis=1)
    values = np.array(df_out['Result'].dropna())
    df_out['z_score']= np.abs(zscore(values))
    df_out.sort_values(by='z_score', ascending=False, inplace=True)
    return df_out

"""Set work directory"""
os.chdir('paste work directory path')
"""File containing gene expression data (gene expression data obtained from CancerRX database)"""
df= pd.read_excel('CancerRX_genes_exp_normalized.xlsx')
"""Data file listing cell types"""
df_cell_types= pd.read_excel('CancerRX_Cell_lines.xlsx')
"""Data set containing drugs' IC50 values from CancerRX database"""
drugs= {'Ulix': 'Ulix_PANCANCER.xlsx',
        'SCH': 'SCH_PANCANCER.xlsx',
        'VX11e': 'VX11e_PANCANCER.xlsx'}
"""Gene sets to analyze"""
gene_sets= {'ERKgenes': 'ERK related genes.xlsx'}
"""Specify which cell type, drug, and gene set to use in analysis"""
cell_type= 'Solid' # Solid or Blood
drug= 'Ulix' # Ulix, SCH or VX11e
gene_set= 'ERKgenes'

df_drug= pd.read_excel(drugs[drug])
df_genes= pd.read_excel(gene_sets[gene_set])
cells= list(set(df['Cell_line_name']) & set(df_drug['Cell line name']) & set(df_cell_types[cell_type])) #selects cell line with gene expression and IC50 data available
df_drug['IC50_exp']= df_drug['IC50'].apply(lambda x: math.exp(x)) #converts IC50 values from ln() to linear values to avoid negative values
df.set_index('Cell_line_name', inplace=True)
df_drug.set_index('Cell line name', inplace=True)
df_comb= pd.concat([df.loc[cells], df_drug.loc[cells]['IC50_exp']], axis=1) #creates DataFrame with combined gene expression and IC50 data
genes= list(set(df_genes['Gene']) & set(df.columns))
Exp=df_comb[genes] #gene expression will be used as an argument
IC50=df_comb['IC50_exp'] #drug IC50 values will be used as predicted value

"""set tolerance parameter which will be used for elastic net"""
tolerance= 0.005
param= elastic_param(Exp, IC50, 10) #estimate best parameters for elastic net 
df_result= perform_elastic(Exp, IC50, param['alpha'], param['l1_ratio'], 100)
"""saves results to Elastic Net subfolder"""
out_dir= os.path.join(os.getcwd(), 'Elastic Net')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
df_result.to_excel(out_dir+'CancerRX_'+cell_type+'_'+drug+'_'+gene_set+'_ElasticNet.xlsx')
