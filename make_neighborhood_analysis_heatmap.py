import pandas as pd
import os
import numpy as np

"""Summarize neighborhood_analysis results in the following steps. This script preprocess the 
data for analysis with Morpheus.
    Mask non-significant neighhood fractions with 0.
    Add metadata for each ROI
    Output a Morpheus ready table for additional analysis with Morpheus.
"""

path = 'N:/HiTS Projects and Data/Personal/Jake/mgh_tma/results'
os.chdir(path)
roi_meta = pd.read_csv('../processed_data/roi_metadata.csv', index_col=0)
site_annotation = pd.read_excel(
    '../processed_data/site_annotation.xlsx', index_col=0)
roi_meta = roi_meta.merge(site_annotation, left_on='Site', right_index=True)
pvals = pd.read_csv('neighborhood_pvalues.csv', index_col=0)
fractions = pd.read_csv('neighborhood_fractions.csv', index_col=0).fillna(0)
# make sure the two tables have the same columns and order
assert((pvals.columns == fractions.columns).all())

pvals.set_index(pvals.group_id + '_' + pvals.index, inplace=True)
fractions.set_index(fractions.group_id + '_' + fractions.index, inplace=True)
# check if indices match
assert((pvals.index == fractions.index).all())

pvals.iloc[:, :-1] = pvals.iloc[:, :-
                                1].apply(lambda x: [True if i <= 0.05 else False for i in x.values])
fractions.iloc[:, :-1] = fractions.iloc[:, :-1] * pvals.iloc[:, :-1]
fractions['cluster'] = ['_'.join(x.split('_')[2:]) for x in fractions.index]
fractions = fractions[~fractions.cluster.isin(
    ['Others', 'CD45_DP', 'CD45_DN'])]
fractions.drop(['Others', 'CD45_DP', 'CD45_DN'], axis=1, inplace=True)

fractions = fractions.reset_index().merge(roi_meta[
    ['group_id', 'patient', 'Site', 'organ']], on='group_id', how='left').set_index('index')
cols = fractions.columns[-5:].tolist() + fractions.columns[:-5].tolist()
fractions = fractions[cols]
fractions.transpose().to_csv('neighborhood_analysis_heatmap.csv')
