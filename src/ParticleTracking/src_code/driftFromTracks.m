function [movieInfo, dfEst] = driftFromTracks(movieInfo,g)
% we use the common motion of adjacent frames to estimate the drifting of
% data

%% change to non_rigid registration
if isfield(g, 'translation_path')  && isfield(g, 'timepts_to_process')
%     file_list = dir(g.translation_path);
%     file_list = file_list(3:end);
    timepts_to_process = g.timepts_to_process(1:end-1);
    grid_size = g.driftInfo.grid_size;
    movieInfo.drift.x_grid = cell(length(timepts_to_process),1);
    movieInfo.drift.y_grid = cell(length(timepts_to_process),1);
    movieInfo.drift.z_grid = cell(length(timepts_to_process),1);
    for ii = 1:length(timepts_to_process)
        load(fullfile(g.translation_path, timepts_to_process(ii)+'.mat'), 'phi_current_vec');
        x_grid = padarray(reshape(phi_current_vec(1:3:end-2), [grid_size grid_size grid_size]), [1 1 1], 'replicate'); 
        y_grid = padarray(reshape(phi_current_vec(2:3:end-1), [grid_size grid_size grid_size]), [1 1 1], 'replicate'); 
        z_grid = padarray(reshape(phi_current_vec(3:3:end), [grid_size grid_size grid_size]), [1 1 1], 'replicate'); 
        movieInfo.drift.x_grid{ii} = x_grid;
        movieInfo.drift.y_grid{ii} = y_grid;
        movieInfo.drift.z_grid{ii} = z_grid;
    end
end
movieInfo.driftInfo = g.driftInfo;  % for conveinent input


%%
movieInfo.xCoord = movieInfo.orgCoord(:,1);
movieInfo.yCoord = movieInfo.orgCoord(:,2);
movieInfo.zCoord = movieInfo.orgCoord(:,3);

t = max(movieInfo.frames);
trackLength = cellfun(@length,movieInfo.tracks);
validTrack = find(trackLength>=g.trackLength4var); % should we?
if length(validTrack) < 10
    dfEst = 0;
    return;
end

% if isfield(g, 'translation') && g.use_translation_as_drift
%     dfEst = g.translation;
% else
%     vMatXYZ = nan(length(validTrack), t-1, 3); % velocity for x,y,z direction
% 
%     for j=1:length(validTrack)
%         trackNum = validTrack(j);
%         curTrack = movieInfo.tracks{trackNum};
%         curFrames = movieInfo.frames(curTrack);
%         adjLinks = find(curFrames(2:end)-curFrames(1:end-1)==1);
%         for k = 1:length(adjLinks)
%             timeStart = curFrames(adjLinks(k));
%             ptStart = curTrack(adjLinks(k));
%             ptEnd = curTrack(adjLinks(k)+1);
%             vMatXYZ(j, timeStart, 1) = movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart);
%             vMatXYZ(j, timeStart, 2) = movieInfo.yCoord(ptEnd)-movieInfo.yCoord(ptStart);
%             vMatXYZ(j, timeStart, 3) = movieInfo.zCoord(ptEnd)-movieInfo.zCoord(ptStart);
%         end
%     end
% 
%     dfEst_perFrame = squeeze(nanmean(vMatXYZ,1));
%     dfEst = zeros(size(dfEst_perFrame,1)+1, size(dfEst_perFrame,2));
%     for i=1:t-1
%         dfEst(i+1,:) = sum(dfEst_perFrame(1:i,:),1);
%     end   
% end
% if ~isfield(g, 'applyDrift2allCoordinate')
%     g.applyDrift2allCoordinate = true; % particle tracking
% end
if g.applyDrift2allCoordinate == true % particle tracking
    % update movieInfo (xCoord, yCoord and zCoord)
    for i=1:t-1
        xyzMotion = dfEst(i+1,:);%sum(dfEst(1:i,:),1);
        uptIdx = movieInfo.frames==(i+1);
        movieInfo.xCoord(uptIdx) = movieInfo.xCoord(uptIdx)-xyzMotion(1);
        movieInfo.yCoord(uptIdx) = movieInfo.yCoord(uptIdx)-xyzMotion(2);
        movieInfo.zCoord(uptIdx) = movieInfo.zCoord(uptIdx)-xyzMotion(3);
    end
    % update vox index
    if isfield(movieInfo, 'vox')
        for i=1:t-1
            xyzMotion = dfEst(i+1,:);%sum(dfEst(1:i,:),1);
            uptIdx = movieInfo.frames==(i+1); % can be optimized
            movieInfo.vox(uptIdx) = cellfun(@(x) x-xyzMotion, movieInfo.vox(uptIdx),...
                'UniformOutput',false);
        end
    end
else
    % update CDist, CDist_i2j, CDist_j2i
    if isfield(movieInfo, 'CDist')
        %vox = movieInfo.vox;
        %nei = movieInfo.nei;
        CDist = cell(numel(movieInfo.CDist), 1);
        CDist_i2j = cell(numel(movieInfo.CDist), 1);
        for i=1:numel(CDist)
            %disp(i);
            cur_Vox = movieInfo.vox{i};
            %cur_nei_vox = movieInfo.vox(movieInfo.nei{i});
            nei_f = movieInfo.frames(movieInfo.nei{i});
            cur_dist = zeros(length(nei_f),1);
            cur_dist_dirWise = zeros(length(nei_f),2);
            for j=1:length(nei_f)
%                 frame_shift = dfEst(nei_f(j),:) - dfEst(movieInfo.frames(i),:);
                frame_shift = getNonRigidDrift(cur_Vox, movieInfo.vox{movieInfo.nei{i}(j)},...
                    movieInfo.frames(i), nei_f(j), movieInfo.drift, g.driftInfo);
                [cur_dist(j), ~, cur_dist_dirWise(j,:)] = ovDistanceRegion(cur_Vox, ...
                    movieInfo.vox{movieInfo.nei{i}(j)}, frame_shift);%
            end
            CDist{i} = cur_dist;
            CDist_i2j{i} = cur_dist_dirWise; % save both i2j and j2i in two columns
        end
    
        CDist_j2i = cell(numel(movieInfo.nei), 1);
        for i=1:numel(movieInfo.nei)
            for j=1:length(movieInfo.nei{i})
                CDist_j2i{movieInfo.nei{i}(j)} = cat(1, CDist_j2i{movieInfo.nei{i}(j)}, ...
                    CDist_i2j{i}(j, 2));
            end
        end
        movieInfo.CDist = CDist;
        movieInfo.CDist_i2j = CDist_i2j;
        movieInfo.CDist_j2i = CDist_j2i;
    end
end
%% change to Mengfan's global registration
% if isfield(g, 'translation_path')
%     global_translation = load(g.translation_path);
%     for i = 1:numel(global_translation.tform)
%         dfEst(i+1,:) = global_translation.tform{i}.T(4, 1:3);
%     end
% end
% movieInfo.drift = dfEst; 
% % neighborhoold should be updated too, but not necessary
% fprintf('finish drifting correction\n');

fprintf('finish drifting correction\n');
