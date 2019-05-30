import pandas as pd
import os

# Process HistoCAT outputs
path = 'N:/HiTS Projects and Data/Personal/Jake/mgh_tma/histocat_output'
os.chdir(path)
done_plate = []
metadata = pd.DataFrame()
if 'expr' in locals():
    del expr
for fn in sorted(os.listdir(path)):
    if 'nucleiMask' not in fn:
        continue
    _metadata = fn.split('_')
    plate = _metadata[1][2:6]
    if plate not in done_plate:
        print(plate, fn)
        if 'expr' in locals():
            expr.to_csv('../processed_data/' +
                        '_'.join([done_plate.pop(0), 'nuclei_log_normed.csv']))
        done_plate.append(plate)
        expr = pd.DataFrame()
    roi = _metadata[2]
    _data_fn = '.'.join([roi, 'csv'])
    _df = pd.read_csv(os.path.join(path, fn, _data_fn), index_col=1)
    _df.index = '_'.join([plate, roi, '']) + _df.index.astype(str)
    _df_expr = _df.loc[:, 'Cell_Marker1':'Cell_Marker44']
    _df_meta = _df.loc[:, 'Area':]
    _df_meta['ROI'] = roi
    _df_meta['Plate'] = plate
    expr = expr.append(_df_expr)
    metadata = metadata.append(_df_meta)
expr.to_csv('../processed_data/' +
            '_'.join([done_plate.pop(0), 'nuclei_log_normed.csv']))
metadata.to_csv('../processed_data/' + 'tma_metadata.csv')
