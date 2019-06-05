import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.spatial.distance import pdist, squareform


def spotwise_clusterdist_plot(plot_data, metric='Euclidean'):
    fig, axes = plt.subplots(1 + plot_data.cluster.nunique() //
                             3, 3, figsize=(16, 16), sharex=True, sharey=True)
    axes = axes.ravel()
    i = 0
    for _group in plot_data.groupby('cluster'):
        group_name, group_df = _group
        group_df = group_df.sort_values(['cluster', 'Site'])

        dist_data = group_df.set_index('Site').iloc[:, 4:]
        dist = squareform(pdist(dist_data, metric=metric))
        site_names = dist_data.index
        dist = pd.DataFrame(dist, index=site_names,
                            columns=site_names).stack()
        dist.index.names = ['Site1', 'Site2']
        dist = dist.reset_index()
        dist = dist[dist.iloc[:, 2] > 0]

        same_spot_dist = dist[dist.Site1 == dist.Site2].iloc[:, 2]
        diff_spot_dist = dist[dist.Site1 != dist.Site2].iloc[:, 2]
        same_spot_dist.name = 'Within_spot_dist'
        diff_spot_dist.name = 'Across_spot_dist'
        sns.kdeplot(same_spot_dist, ax=axes[i])
        sns.kdeplot(diff_spot_dist, ax=axes[i])
        axes[i].set_title(group_name)
        axes[i].set_xlabel(metric + ' distance', fontsize=16)
        i += 1


def check_roi(expr_data, metadata, tma_roi, plotting=False):
    plot_data_idx = metadata[metadata.group_id == tma_roi].index
    if plotting:
        color_col = 'Ecad'
        plot_data = metadata.loc[plot_data_idx]
        expr_data.loc[plot_data_idx, color_col].hist(bins=100)
        x_vector = plot_data.X_position.astype(float).values
        y_vector = plot_data.Y_position.astype(float).values
        x_vector = x_vector - x_vector.min()
        y_vector = y_vector - y_vector.min()
        plot_data[color_col] = expr_data.loc[plot_data_idx, color_col]
        plt.figure(figsize=(9, 9))
        sns.scatterplot(x_vector, y_vector, data=plot_data, hue=color_col,
                        s=2, palette='coolwarm', edgecolor=None, alpha=0.6)
        axis_max = np.max((x_vector.max(), y_vector.max()))
        _ = plt.xlim((0, axis_max))
        _ = plt.ylim((0, axis_max))
        return expr_data.loc[plot_data_idx]
    else:
        return expr_data.loc[plot_data_idx]
