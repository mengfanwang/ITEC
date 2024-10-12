function VarianceMapEstimation(data_path, tmp_path, scale_term, clipping)

tic;
tif_files = dir(fullfile(data_path, '/*.tif'));
if ~exist(tmp_path,'dir')
    mkdir(tmp_path);
end
if ~isfolder(fullfile(tmp_path, 'varianceMap'))
    mkdir(fullfile(tmp_path, 'varianceMap'));
end
for i= 1:numel(tif_files)  
    fprintf('Variance map estimating %d/%d\n', i, numel(tif_files));
    
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    [~, org_name, ~] = fileparts(tif_files(i).name);
    org_im = org_im - clipping;      
    org_im(org_im < 0) = 0;
    vid = 255*org_im/scale_term;
    varMap = cell(3,2);
    [varMap{1,1}, varMap{2,1},varMap{3,1}] = ...
        calVarianceStablizationBY(vid, 0.8, 3);
    vid_stb = sqrt(vid+3/8);
    [varMap{1,2}, varMap{2,2},varMap{3,2}] = ...
        calVarianceStablizationBY(vid_stb, 0.8, 3);
    save(fullfile(tmp_path, 'varianceMap', [org_name '.mat']), 'varMap','-v7.3');
end
fprintf('Variance estimation running time:');
toc

