import pandas as pd
import os
import numpy as np
from sklearn import mixture
from hdbscan import HDBSCAN
from sklearn.cluster import KMeans
from cycifsuite.get_data import read_synapse_file
from sklearn.mixture import GaussianMixture


def iterative_gmm(expr_data, col, cluster_info=None, n_comp=2, iteration='Round_1'):
    """Gating using GMM with 2 or 3 components.

    Parameters
    ========
    expr_data : pd.DataFrame
        log2 transformed data, single cells in rows and marks in columns.
    col : list
        list of columns to run Gaussian mixture models on.
    cluster_info : pd.DataFrame or None
        table containing the cluster information to be updated by the function. 
        if None, a new table will be created with indices the same as the input `expr_data',
        and column as the iteration name.
    iteration : str
        name for the columns in the cluster_info table.

    Returns
    ========
    cluster_info : pd.DataFrame or None
        table containing the cluster information to be updated by the function. 
        if None, a new table will be created with indices the same as the input `expr_data',
        and column as the iteration name.
    """
    gmm = GaussianMixture(n_components=2, tol=1e-9,
                          init_params='random', max_iter=500)
    labels = gmm.fit_predict(expr_data[col])
    marker_cluster_exp = expr_data[col].groupby(labels).median()
    if n_comp == 2:
        descriptor = ['low', 'high']
    elif n_comp == 3:
        descriptor = ['low', 'med', 'high']
    else:
        print('More than 3 components is not implemented, reverting to 3 components')
        descriptor = ['low', 'med', 'high']
        n_comp = 3
    if cluster_info is None:
        cluster_info = pd.DataFrame(index=expr_data.index)
    for cluster, keyword in zip(marker_cluster_exp.sort_values(col).index, descriptor):
        _idx = [i for i, x in enumerate(labels) if x == cluster]
        cluster_info.loc[expr_data.index[_idx],
                         iteration] = col[0] + '_' + keyword
    return cluster_info

if __name__ == '__main__':
    """Hard coded iterative gating, Ecad=>SMA=>CD45, with additional gating on Ki67 and gH2ax.
    """
    path = 'N:/HiTS Projects and Data/Personal/Jake/mgh_tma/processed_data'
    os.chdir(path)
    expr = pd.read_hdf('index_corrected_tma_expr_data.hdf')
    metadata = pd.read_csv('index_corrected_tma_metadata.csv', index_col=0)
    channel_info = pd.read_csv(read_synapse_file(
        'syn18555930'), index_col=0).stack()
    metadata['group_id'] = metadata.Plate + '_' + metadata.ROI.astype(str)
    colnames = channel_info.values
    expr.columns = colnames
    valid_cells = metadata[metadata.labeled_as_lost == 'No'].index
    valid_cols = [x for x in colnames if 'DNA' not in x]
    expr = expr.reindex(valid_cells)[valid_cols]
    metadata = metadata.loc[valid_cells]

    cluster_info = pd.DataFrame()
    for plate_id in ['TMA1', 'TMA2', 'TMA3', 'TMA4']:
        plate_meta = metadata[(metadata.Plate == plate_id) & (
            metadata.labeled_as_lost == 'No')][['ROI']]
        patient_meta = pd.read_excel(
            '../CMTMA_Breast_CDK4 autopsies_DEIDENTIFIED_JRL_20181119.xlsx', sheet_name=plate_id, index_col=0, dtype=str)
        plate_meta.ROI = plate_meta.ROI.astype('int64')
        plate_meta = plate_meta.merge(
            patient_meta.iloc[:, :4], left_on='ROI', right_index=True, how='left')
        plate_meta['patient'] = [x.split('-')[0] for x in plate_meta.RAN_UNI]
        for group in plate_meta.groupby('patient'):
            patient, patient_df = group
            patient_expr = expr.loc[patient_df.index]
            print('processing patient: {}'.format(patient))
            patient_plate_meta = plate_meta[plate_meta.patient == patient]
            # patient_umap = dimension_reduction(patient_expr, umap_n_neighbor_modifier=1000)

            # Round 1
            clustered = iterative_gmm(
                patient_expr, ['Ecad', 'CK8-FITC'],  iteration='Ecad')

            # Round 2-a
            _current_round_data = patient_expr.loc[
                clustered[clustered.Ecad == 'Ecad_low'].index]
            round_name = 'Round_2a'
            clustered = iterative_gmm(_current_round_data, [
                                      'aSMA'], cluster_info=clustered, iteration=round_name)

            # Round 2-b
            _current_round_data = patient_expr.loc[
                clustered[clustered.Ecad == 'Ecad_high'].index]
            round_name = 'Round_2b'
            clustered = iterative_gmm(_current_round_data, [
                                      'gH2ax-PE'], cluster_info=clustered, iteration=round_name)

            # Round 3
            _current_round_data = patient_expr.loc[
                clustered[clustered.Round_2a == 'aSMA_low'].index]
            round_name = 'Round_3'
            clustered = iterative_gmm(_current_round_data, [
                                      'CD45-PE'], cluster_info=clustered, iteration=round_name)

            # Round 4
            _current_round_data = patient_expr.loc[
                clustered[clustered.Round_3 == 'CD45-PE_high'].index]
            clustered = iterative_gmm(_current_round_data, [
                                      'CD4'], cluster_info=clustered, iteration='CD4')
            clustered = iterative_gmm(_current_round_data, [
                                      'CD8a'], cluster_info=clustered, iteration='CD8')

            # Ki67 overall
            clustered = iterative_gmm(
                patient_expr, ['Ki67-570'], cluster_info=clustered, iteration='Ki67')

            # Finalize cluster names
            clustered.fillna('', inplace=True)

            # Add patient info
            clustered['patient'] = patient
            cluster_info = cluster_info.append(clustered)

    cluster_info['cluster_name'] = cluster_info.loc[
        :, 'Ecad':'Ki67'].apply(lambda x: '|'.join(x.values), axis=1)
    new_cluster_names = ['Epi|gH2_high|Ki67_high', 'Epi|gH2_high|Ki67_low', 'Epi|gH2_low|Ki67_high', 'Epi|gH2_low|Ki67_low',
                         'Stromal|Ki67_high', 'Stromal|Ki67_low',
                         'CD45_DP', 'CD45_DP', 'CD45_CD4', 'CD45_CD4',
                         'CD45_CD8', 'CD45_CD8', 'CD45_DN', 'CD45_DN', 'Others', 'Others']
    for old_cluster_name, new_name in zip(sorted(cluster_info.cluster_name.unique()), new_cluster_names):
        cluster_info.cluster_name.replace(
            old_cluster_name, new_name, inplace=True)
    # Update metadata
    metadata.loc[cluster_info.index, 'cluster'] = cluster_info[
        'cluster_name'].values
    metadata.to_csv('../results/clustered_metadata.csv')
