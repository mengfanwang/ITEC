function orgIm3d = imdisplayWithTrack3D(orgIm3d, movieInfo, timePt, p)
% generate a 3D colorful data with ROI overlay
particleSize = p.particleSize;
%particleCl = p.particleCl;
stoppedParticleCl = p.stoppedParticleCl;
lineWidth = p.lineWidth;
RegLineCl = p.RegLineCl;
stoppeLineCl = p.stoppeLineCl;
clMap = p.cmap;
for i=1:numel(movieInfo.tracks)
    particleCl = clMap(i,:);
    fms = movieInfo.frames(movieInfo.tracks{i});
    % continuous path
    node = find(fms==timePt);% at most one
    if ~isempty(node)
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
                if j==node
                    orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, particleCl);
                end
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
                    if j==node
                        orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, stoppedParticleCl);
                    end
                    orgIm3d = draw3Dline(orgIm3d, hPt, tPt, lineWidth, stoppeLineCl);
                end
            end
        end
    end
end


end
