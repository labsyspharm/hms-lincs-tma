import pandas as pd
import os
from cycifsuite.detect_lost_cells import *
from cycifsuite.common_apis import channel_histograms

path = 'N:/HiTS Projects and Data/Personal/Jake/mgh_tma/processed_data'
os.chdir(path)
plates = ['TMA' + str(x) + '_nuclei_log_normed.csv' for x in range(1, 5)]
metadata = pd.read_csv('tma_metadata.csv', index_col=0)
dna_channels = np.arange(0, 44, 4)[1:]
total_loss = []
for plate in plates:
    fn = plate
    expr = pd.read_csv(fn, index_col=0)
    expr_qc = expr.iloc[:, dna_channels]
    expr_qc = 2**expr_qc
    x, y, t1 = ROC_lostcells(expr_qc, n_cycles=10, cutoff_min=0, cutoff_max=1, steps=50,
                             segmentation_cycle=1, filtering_method='cycle_diff', left_stepping=False,
                             figname='../results/' + plate.replace('csv', 'png'))
    _, lc_cv_down = get_lost_cells(expr_qc, t1, 10, 'cycle_diff')
    total_loss += list(set(lc_cv_down))
metadata = update_metadata(metadata, total_loss)

# Nuclei size thresholding
threshold = 9
metadata.loc[metadata.Area < threshold, 'labeled_as_lost'] = 'Yes'
metadata.to_csv('tma_metadata.csv')
