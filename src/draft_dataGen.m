clc;clear;close all;
% the sample data is cropped from 21-01-11, 130-149tps, left uppper corner
% (600,320,130), sz = [800 800 25].

org_path = '/work/Mengfan/Embryo/22-01-11/sameViewFusion_10';
tgt_path = '../data/tiff';

sz = [800 800 25];
for ii = 0:19
    ii
    im = tifread(fullfile(org_path, [sprintf('%05d', ii+130) '.tif']));
    im = im(601:600+sz(1), 321:320+sz(2), 131:130+sz(3));
    tifwrite(uint16(im), fullfile(tgt_path, sprintf('%05d', ii)));
end


%% convert to hdf5file
timepts_to_process = generate_tps_str(0:19);
tif2bdv(tgt_path, '../data/hdf5/embryo_data_h5', timepts_to_process, [], []);