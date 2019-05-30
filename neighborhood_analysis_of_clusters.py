import pandas as pd
import numpy as np
import os


def neighborhood_analysis(metadata, grouping_col='group_id', cluster_col='cluster', **kwargs):
    """Overall script for neighborhood analysis.
    Each spot is analyzed separately. Neighbors of each cluster within were evaluated individually
    and sequentially. Neighbor cells of each single cell within a cluster are pooled and considered neighbors
    of that particular cluster. Pvalues estimated by permutation analysis and the raw fractions are exported.

    Parameters
    ========
    metadata : pd.DataFrame
        the table containing cluster in formation and spot information. The metadata table must have `neighbour`
        columns from histoCAT analysis.
    (grouping, cluster)_col : str
        column names in the metadata table containing the spot information and the cluster information.
    **kwargs : additional arguments in 'permutation_neighborhood' function, where 'num_permutations' controls
        numbers of permutations.

    Returns
    ========
    pval_report : pd.DataFrame
        a table containing all estimated pvalues of getting the observed neighbor fraction by chance. 
    fraction_report : pd.DataFrame
        a table containing all observed neighborhood fractions. 
    """

    pval_report = pd.DataFrame()
    fraction_report = pd.DataFrame()
    for _group in metadata.groupby(grouping_col):
        group_name, group_df = _group
        group_pvals = pd.DataFrame()
        group_fractions = pd.DataFrame()
        for _cluster_group in group_df.groupby(cluster_col):
            cluster_name, cluster_df = _cluster_group
            cluster_neighbor_pvals, _neighbor_fractions = permutation_neighborhood(
                cluster_df, group_df, **kwargs)
            cluster_neighbor_pvals = pd.DataFrame(
                cluster_neighbor_pvals, columns=[cluster_name])
            group_pvals = pd.concat(
                [group_pvals, cluster_neighbor_pvals], axis=1, sort=False)
            _neighbor_fractions.name = cluster_name
            group_fractions = pd.concat(
                [group_fractions, _neighbor_fractions], axis=1, sort=False)
        # record p-values
        group_pvals[grouping_col] = group_name
        group_pvals.fillna(1, inplace=True)
        pval_report = pval_report.append(group_pvals, sort=False)
        # record detailed neighbor fractions
        group_fractions.fillna(0, inplace=True)
        group_fractions[grouping_col] = group_name
        fraction_report = fraction_report.append(group_fractions, sort=False)
    return pval_report, fraction_report


def get_neighbors(neighbour_metadata, exclude_self=True):
    neighbour_cols = [
        x for x in neighbour_metadata.columns if 'neighbour' in x]
    spot_prefix = '_'.join(neighbour_metadata.index[0].split('_')[:2])
    cell_ids = neighbour_metadata[neighbour_cols].values.flatten()
    cell_ids = cell_ids[(~np.isnan(cell_ids)) & (cell_ids != 0)]
    cell_ids = np.unique(cell_ids).astype(int)
    neighbors_idx = [spot_prefix + '_' + str(x) for x in cell_ids]
    if exclude_self:
        neighbour_cols = [
            x for x in neighbors_idx if x not in neighbour_metadata.index]
    return neighbors_idx


def annotate_neighbors(neighbour_idx, metadata):
    neighbors = metadata.reindex(neighbour_idx).dropna(
        how='all').cluster.value_counts()
    total_neighbors = neighbors.sum()
    neighbour_fractions = neighbors / total_neighbors
    return neighbour_fractions.sort_index()


def neighbor_across_spots(metadata, target_cluster, spot_col='group_id', cluster_col='cluster'):
    """Get neighbor information for the target cluster across all the spots available in the metadata sheet.
    """
    all_spots = metadata[spot_col].unique()
    neighbor_summary = pd.DataFrame()
    for _group in metadata.groupby(spot_col):
        spot_name, spot_metadata = _group
        spot_metadata = spot_metadata[spot_metadata.cluster == target_cluster]
        if spot_metadata.shape[0] == 0:
            continue
        cluster_neighbors = get_neighbors(spot_metadata)
        annotated_neighbors = annotate_neighbors(cluster_neighbors, metadata)
        annotated_neighbors.name = spot_name
        annotated_neighbors = pd.DataFrame(annotated_neighbors).transpose()
        neighbor_summary = neighbor_summary.append(annotated_neighbors)
    return neighbor_summary.fillna(0)


def permutation_neighborhood(target_metadata, spot_metadata, num_permutations=1000, verbose=False):
    """Permutate cluster labels in the spot_metadata table to get a null distribution of observing the neighbor by cluster profile by chance. 
    """
    neighbors = get_neighbors(target_metadata, exclude_self=False)
    true_neighbor_fractions = annotate_neighbors(neighbors, spot_metadata)
    pval = pd.Series(0, index=true_neighbor_fractions.index)
    for i in range(num_permutations):
        permuted_cluster = np.random.permutation(spot_metadata.cluster.values)
        permuted_data = spot_metadata.copy()
        permuted_data.loc[:, 'cluster'] = permuted_cluster
        permutated_neighbor_fractions = annotate_neighbors(
            neighbors, permuted_data)
        permutated_neighbor_fractions = permutated_neighbor_fractions.reindex(
            true_neighbor_fractions.index, fill_value=0)
        pval += (true_neighbor_fractions <= permutated_neighbor_fractions)
        if (i % 100 == 0) & (verbose):
            print('Iteration: {}'.format(str(i)))
            print(permutated_neighbor_fractions)
    return pval / num_permutations, true_neighbor_fractions

if __name__ == '__main__':
    """Hard coded neighborhood analysis with default parameters.
    """
    path = 'N:/HiTS Projects and Data/Personal/Jake/mgh_tma/processed_data'
    os.chdir(path)
    metadata = pd.read_csv('clustered_metadata.csv', index_col=0)
    metadata['group_id'] = metadata.Plate + '_' + metadata.ROI.astype(str)
    pvals, fractions = neighborhood_analysis(metadata)
    pvals.to_csv('../results/neighborhood_pvalues.csv')
    fractions.to_csv('../results/neighborhood_fractions.csv')
