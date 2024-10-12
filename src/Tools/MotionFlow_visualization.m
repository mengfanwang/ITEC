clc;clear;close all;

data_folder = '../../data/tiff';
tmp_path = '../../tmp';
reg_ind = 0;

data1 = tifread(fullfile(data_folder, ['0000' num2str(reg_ind) '.tif']));
fore2 = load(fullfile(tmp_path, 'InstanceSeg_res', ['0000' num2str(reg_ind) '.mat']));
fore2 = fore2.refine_res>0;
load(fullfile(tmp_path, 'MotionFlow', ['0000' num2str(reg_ind) '.mat']));
data1_backup = data1;
fore2_backup = fore2;

%% plot
zz = 10;
data1 = data1_backup(:,:,zz)/2;
fore2 = fore2_backup(:,:,zz);

phi_current_vec = reshape(phi_current_vec,3,[]);
phi_current = imresize4d(phi_current_vec, [800 800 32], 'linear');
ux = phi_current(:,:,zz,1);
uy = phi_current(:,:,zz,2);

ux(~fore2) = 0;
uy(~fore2) = 0;
uz(~fore2) = 0;


[x_size, y_size, z_size] = size(ux);
% Enhance the quiver plot visually by downsizing vectors  
%   -f : downsizing factor
v = sqrt(ux.^2 + uy.^2);
% x = zeros(x_size, y_size); y = zeros(x_size, y_size);
% [~, v_maxind] = max(v,[],3);
% for ii = 1:size(v)
%     for jj = 1:size(v)
%         x(ii,jj) = ux(ii,jj,v_maxind(ii,jj));
%         y(ii,jj) = uy(ii,jj,v_maxind(ii,jj));
%     end
% end
f=8;
color = colormap;
v = v(1:f:x_size,1:f:y_size);
color_quiver = ceil(sqrt(max(v,[],3))*255/sqrt(max(v(:))+1e-6)) + 1; % sqrt: change vector color
x = ux(1:f:x_size,1:f:y_size);
y = uy(1:f:x_size,1:f:y_size);
[X,Y]=meshgrid(1:size(x,2),1:size(x,1));
scale = f;
h = imshow(data1/255);hold on;
% h = imshow(imresize(max((data1-data2).^2*2/3,[],3),4));hold on;
% h = imshow(zeros(size(imresize(max(data1,[],3)/2,scale)))); hold on;
for ii = 1:256
    ind = color_quiver == ii;
    quiver(X(ind)*scale,Y(ind)*scale,-y(ind),-x(ind),0,'Color',color(ii,:));
end
% camroll(-90);

function phi_current = imresize4d(phi_current, scale, method)
    % transform to 4d data if table
    % resize 4d data
    % gpu not supported
    flag = 0;
    if isgpuarray(phi_current)
        isgpu = 1;
        phi_current = gather(phi_current);
    else
        isgpu = 0;
    end
    if ismatrix(phi_current)       % if table
        % change to [3 []]
        if size(phi_current,1) == 3
            flag = 1;
        elseif size(phi_current, 2) == 3
            flag = 1;
            phi_current = phi_current';
        end
        if flag
            batch_dim = round(size(phi_current,2)^(1/3));
            if batch_dim^3 ~= size(phi_current,2)
                error('Incorrect number of elements.');
            end
            if length(scale) == 1
                scale = [batch_dim batch_dim batch_dim] * scale;
            end
            phi_current_temp = phi_current';
            phi_current = zeros([scale 3]);
            if size(phi_current_temp,1) == 1
                phi_current(:,:,:,1) = phi_current_temp(:,1);
                phi_current(:,:,:,2) = phi_current_temp(:,2);
                phi_current(:,:,:,3) = phi_current_temp(:,3);
            else
                phi_current(:,:,:,1) = imresize3(reshape(phi_current_temp(:,1), batch_dim, batch_dim, batch_dim), scale, method);
                phi_current(:,:,:,2) = imresize3(reshape(phi_current_temp(:,2), batch_dim, batch_dim, batch_dim), scale, method);
                phi_current(:,:,:,3) = imresize3(reshape(phi_current_temp(:,3), batch_dim, batch_dim, batch_dim), scale, method);
            end
        end
    elseif  ndims(phi_current) == 4
        if size(phi_current,4) == 3
            flag = 1;
            if length(scale) == 1
                scale = size(phi_current(:,:,:,1)) * scale;
            end
            phi_current_temp = phi_current;
            phi_current = zeros([scale 3]);
            phi_current(:,:,:,1) = imresize3(phi_current_temp(:,:,:,1), scale, method);
            phi_current(:,:,:,2) = imresize3(phi_current_temp(:,:,:,2), scale, method);
            phi_current(:,:,:,3) = imresize3(phi_current_temp(:,:,:,3), scale, method);
        end
    end
    if isgpu
        phi_current = gpuArray(phi_current);
    end
    if flag == 0
        warning('Incorrect input. No operations.');
    end
    
end