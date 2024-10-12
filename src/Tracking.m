function Tracking(data_folder, tmp_path, save_folder, tps_to_process, q)
%% INPUT:
% tps_to_process: specify time points for processing. File names must be
%                 five digits. Set to empty to process all files.

%% define paths
timepts_to_process = generate_tps_str(tps_to_process);
cell_seg_res_folder = tmp_path;
refine_res_folder = fullfile(tmp_path, 'InstanceSeg_res');
if ~exist(save_folder,'dir')
    mkdir(save_folder);
end

%% data preparation
% First downsample the detection results if any
sc_f = q.scaling; % we resize the data to [h/sc_f, w/sc_f, z, t]
st_loc = q.stLoc;
sz_crop = q.szCrop;

% get file name
tif_files = dir(fullfile(data_folder, '*.tif'));
if ~isempty(timepts_to_process)
    tif_files(numel(timepts_to_process)+1:end) = [];
    for f = 1:numel(timepts_to_process)
        tif_files(f).name = timepts_to_process(f) + '.tif';
    end
end

% read synquant res
file_num = numel(tif_files);
[refine_res, threshold_res] = matfiles2cell_scaling(refine_res_folder, ...
    'synQuant_refine_res', timepts_to_process, sc_f, st_loc, sz_crop);
if file_num ~= numel(refine_res)
    warning("The number of segmentation results does not equal to number of files!");
end

% read data
embryo_vid = cell(numel(tif_files), 1); % original data
for i=1:numel(tif_files)
    embryo_vid_temp{1} = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    embryo_vid_temp{1} = 255*embryo_vid_temp{1}./q.scaleTerm;
    [~, embryo_vid_temp, ~, ~, ~, ~] = data_scaling(sc_f, st_loc, ...
    sz_crop, {}, embryo_vid_temp, {}, {}, {}, {});
    embryo_vid{i} = embryo_vid_temp{1};
end
clear embryo_vid_temp

% read variance map
[varMaps, ~] = matfiles2cell_scaling(fullfile(cell_seg_res_folder, 'varianceMap'), ...
    'varianceMap', timepts_to_process, sc_f, st_loc, sz_crop);

% read pricinple curvature result
[eig_res_2d, eig_res_3d] = matfiles2cell_scaling(fullfile(cell_seg_res_folder, 'priCvt'), ...
    'priCvt', timepts_to_process, sc_f, st_loc, sz_crop);
eigMaps = cell(file_num,1);
for i=1:file_num
    eigMaps{i} = cell(2,1);
    eigMaps{i}{1} = eig_res_2d{i};
    eigMaps{i}{2} = eig_res_3d{i};
    eig_res_2d{i} = [];
    eig_res_3d{i} = [];
end

%% parameter setting
g = graphPara_cell(sum(cellfun(@(x) max(x(:)), refine_res)));%1
g.timepts_to_process = timepts_to_process;
g.translation_path = fullfile(tmp_path, 'MotionFlow');
g.driftInfo.grid_size = q.grid_size;    % vector numbers each dim, currently cube only
g.driftInfo.batch_size = q.batch_size;

if g.applyDrift2allCoordinate
    error('Non-rigid version doesn''t accept the option.');
end
[y_batch, x_batch, z_batch] = meshgrid(0:g.driftInfo.grid_size+1);
if ~isempty(st_loc)  % adjust if crop
    y_batch = y_batch*g.driftInfo.batch_size(2) + 0.5 - g.driftInfo.batch_size(2)/2 - st_loc(2)/sc_f;
    x_batch = x_batch*g.driftInfo.batch_size(1) + 0.5 - g.driftInfo.batch_size(1)/2 - st_loc(1)/sc_f;
    z_batch = z_batch*g.driftInfo.batch_size(3) + 0.5 - g.driftInfo.batch_size(3)/2 - st_loc(3);
else
    y_batch = y_batch*g.driftInfo.batch_size(2) + 0.5 - g.driftInfo.batch_size(2)/2;
    x_batch = x_batch*g.driftInfo.batch_size(1) + 0.5 - g.driftInfo.batch_size(1)/2;
    z_batch = z_batch*g.driftInfo.batch_size(3) + 0.5 - g.driftInfo.batch_size(3)/2;
end
g.driftInfo.y_batch = y_batch;
g.driftInfo.x_batch = x_batch;
g.driftInfo.z_batch = z_batch;  % vector locations in original image with one padding

%% start tracking    
q.save_folder = save_folder;
diary(fullfile(save_folder, 'log'));

[movieInfo, movieInfoAll, out_refine_res, refine_resAll,...
    threshold_res, threshold_resAll] = .... 
    mcfTracking_cell(refine_res, embryo_vid, threshold_res, ...
    varMaps, eigMaps, g, q);

%% save results
if q.saveInterMediateRes
    save(fullfile(save_folder, 'movieInfo.mat'), 'movieInfo',...
        'movieInfoAll','-v7.3');
    save(fullfile(save_folder, 'refine_res.mat'), 'out_refine_res',...
        'refine_resAll','-v7.3');
    save(fullfile(save_folder, 'threshold_res.mat'), 'threshold_res',...
        'threshold_resAll','-v7.3');
else
    save(fullfile(save_folder, 'movieInfo.mat'), 'movieInfo', '-v7.3');
    save(fullfile(save_folder, 'refine_res.mat'), 'out_refine_res','-v7.3');
    save(fullfile(save_folder, 'threshold_res.mat'), 'threshold_res','-v7.3');
end


%% save results to mastodon
mastodon_dir = fullfile(save_folder, 'mastodon');
if ~exist(mastodon_dir)
    mkdir(mastodon_dir);
end
addpath('TGMM_wrapper/');
mat2tgmm(movieInfo, fullfile(mastodon_dir, 'tgmm_format'));
tif2bdv(data_folder, fullfile(mastodon_dir, 'embryo_data_h5'), timepts_to_process, st_loc, sz_crop);
