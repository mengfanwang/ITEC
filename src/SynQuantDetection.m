function SynQuantDetection(data_path, tmp_path, q)

tif_files = dir(fullfile(data_path, '/*.tif'));
if ~exist(tmp_path,'dir')
    mkdir(tmp_path);
end
if q.minIntensity < q.clipping
    error('minIntensity must larger than background intensity!');
end
q.minIntensity = q.minIntensity - q.clipping;

%% synQuant
tic;
% add synQuant java path
Pij = fullfile('./src_code_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile('./src_code_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);
p0 = fullfile('./src_code_synquant/SynQuantVid_v1.2.5.1.jar');
javaaddpath(p0);

if ~isfolder(fullfile(tmp_path, 'synQuant_res'))
    mkdir(fullfile(tmp_path, 'synQuant_res'));
end
if ~isfolder(fullfile(tmp_path, 'synQuant_res_tif')) && q.visualization
    mkdir(fullfile(tmp_path, 'synQuant_res_tif'));
end
ds_scale = q.downSampleScale; % down sample scale
for i=1:numel(tif_files)
    fprintf('SynQuant processing %d/%d\n', i, numel(tif_files));
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    [~, org_name, ~] = fileparts(tif_files(i).name);
    [h, w, slices] = size(org_im);
    org_im = imresize3(org_im,round([h/ds_scale w/ds_scale slices]));
    org_im = org_im - q.clipping;  
    org_im(org_im < 0) = 0;
   
    sm_im = imgaussfilt3(org_im,q.sigma);
    [zMap, synId, fMap] = m_Synquant4Embryo_Paramater(sm_im, q);
        
    z_mat = single(zMap);
    id_mat = uint16(synId);
    fMaps = fMap;

    z_mat = imresize3(z_mat,[h w slices],'nearest');
    id_mat = imresize3(id_mat,[h w slices],'nearest');
    fMaps = imresize3(fMaps,[h w slices],'nearest');
    
    save(fullfile(tmp_path, 'synQuant_res', [org_name '.mat']), 'z_mat', 'id_mat','fMaps','-v7.3');
    if q.visualization
        org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
        org_im = org_im - q.clipping;  
        org_im(org_im < 0) = 0;
        labelwrite(uint8(org_im/2), id_mat, fullfile(tmp_path, 'synQuant_res_tif', org_name));
    end
end

% remove java path
javarmpath(p0);
javarmpath(p1);
javarmpath(Pij);
fprintf('Synquant part running time:') 
toc