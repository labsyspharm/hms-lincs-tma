overall_path = 'Z:/sorger/data/RareCyte/089-CM_MGHAUTOPSYTMA-2018APR/ashlar-20190111/CMTMA41/dearray_CMTMA41/';
for n=1:42
    tif_fn = strcat(num2str(n),'.tif');
    mask_dir = strcat(overall_path,'output/',num2str(n),'/');
    % for mask_type={'cellMask.tif', 'cytoplasmMask.tif','cytoplasmOutlines.tif','nucleiMask.tif','nucleiOutlines.tif'}
    for mask_type={'nucleiMask.tif'}
        mask_fn = strcat(num2str(n),'_',mask_type);
        mask_fn = mask_fn{1,1};
        Headless_histoCAT_loading(overall_path,tif_fn, mask_dir,mask_fn, 'd:/histocat/tma_histocat_channel.csv');
    end
end