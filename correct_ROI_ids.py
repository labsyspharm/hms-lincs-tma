# Correct ROI IDs

import pandas as pd
import os

path = 'N:/HiTS Projects and Data/Personal/Jake/mgh_tma/histocat_output'
os.chdir(path)
ashlar_meta = pd.read_excel(
    '../processed_data/tma_ROI ashlar mapping.xlsx', sheet_name=0, skiprows=2, index_col=0)
metadata = pd.read_csv('../processed_data/tma_metadata.csv', index_col=0)
new_metadata = pd.DataFrame()
expr = pd.DataFrame()
for i in range(1, 5):
    print('processing plate {}'.format(str(i)))
    fn = 'TMA0_nuclei_log_normed.csv'.replace('0', str(i))
    _expr = pd.read_csv('../processed_data/' + fn, index_col=0)
    _expr_meta = metadata.loc[_expr.index]
    _expr_meta['real_ROI'] = ashlar_meta.reset_index().set_index(
        'TMA{}'.format(str(i))).loc[_expr_meta.ROI.values, 'Ashlar_ROI'].values
    new_indx = _expr_meta.apply(lambda x: '_'.join(
        [x.name.split('_')[0], str(x.real_ROI), x.name.split('_')[2]]), axis=1).values
    _expr_meta.index = new_indx
    _expr.index = new_indx
    expr = expr.append(_expr)
    new_metadata = new_metadata.append(_expr_meta)
expr.to_hdf('../processed_data/index_corrected_tma_expr_data.hdf',
            key='corrected')
new_metadata.ROI = new_metadata.real_ROI
new_metadata.drop('real_ROI', axis=1, inplace=True)
new_metadata.to_csv('../processed_data/index_corrected_tma_metadata.csv')
