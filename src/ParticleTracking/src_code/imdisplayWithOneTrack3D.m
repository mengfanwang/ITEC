function orgIm3d = imdisplayWithOneTrack3D(orgIm3d, movieInfo, trackIdx, timePt, p, leftUp)
% generate a 3D colorful data with only one Track
% can be merged to imdisplayWithTrack3D
if nargin < 6
    leftUp = [0 0 0];
end
maxInt = p.maxInt;%max(orgIm3d(:));
particleSize = p.particleSize;
particleCl = p.particleCl ;
stoppedParticleCl = p.stoppedParticleCl;
lineWidth = p.lineWidth;
RegLineCl = p.RegLineCl;
stoppeLineCl = p.stoppeLineCl;
[h,w,~,z] = size(orgIm3d);
for i=trackIdx 
    cur_frames = movieInfo.frames(movieInfo.tracks{i});
    uni_frs = unique(cur_frames);
    for j=1:length(uni_frs)
        nodes_cur_fr = movieInfo.tracks{i}(cur_frames==uni_frs(j));
        if length(nodes_cur_fr)>1
            movieInfo.vox{nodes_cur_fr(1)} = cat(1, movieInfo.vox{nodes_cur_fr});
            movieInfo.voxIdx{nodes_cur_fr(1)} = cat(1, movieInfo.voxIdx{nodes_cur_fr});
            movieInfo.tracks{i}(cur_frames==uni_frs(j)) = nodes_cur_fr(1);
            movieInfo.yCoord(nodes_cur_fr(1)) = mean(movieInfo.yCoord(nodes_cur_fr));
            movieInfo.xCoord(nodes_cur_fr(1)) = mean(movieInfo.xCoord(nodes_cur_fr));
            movieInfo.zCoord(nodes_cur_fr(1)) = mean(movieInfo.zCoord(nodes_cur_fr));
        end
    end
    movieInfo.tracks{i} = unique(movieInfo.tracks{i});
end
for i=trackIdx % this is the only difference from function imdisplayWithTrack3D
    fms = movieInfo.frames(movieInfo.tracks{i});
    
    % continuous path
    node = find(fms==timePt);% at most one
    if ~isempty(node)
        if isfield(movieInfo,'voxIdx')
            cur_id = movieInfo.tracks{i}(node);
            voxes = movieInfo.vox{cur_id};
            if isfield(movieInfo,'drift')
                voxes = round(voxes + movieInfo.drift(timePt,:));
            end
            tmp = voxes(:,1);
            voxes(:,1) = voxes(:,2);
            voxes(:,2) = tmp;
            voxes = voxes - leftUp;
            in_val = voxes(:,1)<1 | voxes(:,1)>h | voxes(:,2)<1 | voxes(:,2)>w | voxes(:,3)<1 | voxes(:,3)>z;
            voxIdx = sub2ind([h,w,z], voxes(~in_val,1),voxes(~in_val,2),voxes(~in_val,3));
            %voxIdx = movieInfo.vox{cur_id};
            r = squeeze(orgIm3d(:,:,1,:));
            g = squeeze(orgIm3d(:,:,2,:));
            b = squeeze(orgIm3d(:,:,3,:));
            center_idx = find(b>r); % find the blue center labels
            %r(voxIdx) = 0;
            g(voxIdx) = 0;
            center_vals = b(center_idx);
            b(voxIdx) = 0;%particleCl(3);
            b(center_idx) = center_vals;
            for j=1:size(r,3)
                orgIm3d(:,:,1,j) = r(:,:,j);
                orgIm3d(:,:,2,j) = g(:,:,j);
                orgIm3d(:,:,3,j) = b(:,:,j);
            end
        end
        if length(movieInfo.tracks{i})==1
            pt = [movieInfo.yCoord(movieInfo.tracks{i}), movieInfo.xCoord(movieInfo.tracks{i}), movieInfo.zCoord(movieInfo.tracks{i})];
            orgIm3d = draw3Dparticle(orgIm3d,  pt, particleSize, particleCl);
        elseif timePt==1 || node==1
            head = movieInfo.tracks{i}(1);
            hPt = [movieInfo.yCoord(head),movieInfo.xCoord(head),movieInfo.zCoord(head)];
            orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, particleCl);
        else
            for j = node:-1:2
                tail = movieInfo.tracks{i}(j-1);
                head = movieInfo.tracks{i}(j);
                tPt = [movieInfo.yCoord(tail),movieInfo.xCoord(tail), movieInfo.zCoord(tail)];
                hPt = [movieInfo.yCoord(head),movieInfo.xCoord(head),movieInfo.zCoord(head)];
                orgIm3d = draw3Dparticle(orgIm3d,  tPt, particleSize, particleCl);
                orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, particleCl);
                orgIm3d = draw3Dline(orgIm3d, hPt, tPt, lineWidth, RegLineCl);
            end
        end
    else
        % if the path is jumping this frame
        if sum(fms>timePt) && sum(fms<timePt)
            node = find(fms<timePt);
            node = node(end);
            if node==1
                head = movieInfo.tracks{i}(1);
                hPt = [movieInfo.yCoord(head),movieInfo.xCoord(head),movieInfo.zCoord(head)];
                orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, stoppedParticleCl);
            else
                for j = node:-1:2
                    tail = movieInfo.tracks{i}(j-1);
                    head = movieInfo.tracks{i}(j);
                    tPt = [movieInfo.yCoord(tail),movieInfo.xCoord(tail), movieInfo.zCoord(tail)];
                    hPt = [movieInfo.yCoord(head),movieInfo.xCoord(head),movieInfo.zCoord(head)];
                    orgIm3d = draw3Dparticle(orgIm3d,  tPt, particleSize, stoppedParticleCl);
                    orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, stoppedParticleCl);
                    orgIm3d = draw3Dline(orgIm3d, hPt, tPt, lineWidth, stoppeLineCl);
                end
            end
        end
    end
end

end