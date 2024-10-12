function [v2v1Var, v0Var, s] = v2v1VarCal(movieInfo, g)
% get the variance between (observed velcoity - estimated velocity)
% we have multiple conditions related to k, preNum, postNum

xCoord = movieInfo.xCoord;
yCoord = movieInfo.yCoord;
zCoord = movieInfo.zCoord;
numTrajectories = numel(movieInfo.tracks);
maxEdge = g.validPre+g.validPost;
vCell = cell(maxEdge+1, maxEdge+1);
v2v1Var= nan((maxEdge+1)^2, 10);% preNum, postNum and jumpNum
v2v1Cnt = 1;

tmpG = g;
for preNum = 0:maxEdge
    %disp(preNum)
    for postNum = 0:maxEdge
%         if ~((preNum==1 && postNum==1) ...
%                 || (preNum==g.validPre && postNum==g.validPost))
%             v2v1Cnt = v2v1Cnt+1;
%             continue;
%         end
        if preNum==0 && postNum==0
            v2v1Cnt = v2v1Cnt+1;
            continue;
        end
        tmpG.validPre = preNum;
        tmpG.validPost = postNum;
        tmpG.trackLength4var = max(g.trackLength4var, preNum+postNum+2);
        % real velocity, est v, preNum,postNum, jumpNum, observation variance
        vCell{preNum+1, postNum+1} = nan(1e6, 12);
        cnt = 1;
        for i=1:numTrajectories
            curTrack = movieInfo.tracks{i};
            if length(curTrack)<tmpG.trackLength4var
                continue;
            end
            curTps = movieInfo.frames(curTrack);
            xVec = xCoord(curTrack);
            yVec = yCoord(curTrack);
            zVec = zCoord(curTrack);
            Vec = [xVec, yVec, zVec];
            %curVar = var((Vec(2:end,:)-Vec(1:end-1,:))./(curTps(2:end)-curTps(1:end-1)));
            % estimated velocity
            for j=1:length(curTrack)-1
                if curTps(j+1)-curTps(j)~=1
                    continue;
                end
%                 if curTps(j+1)>30
%                     continue;
%                 end
                [tmpPtsPre, tmpPtsPost] = determinePrePostPts(j, j+1, Vec, tmpG);
                if length(tmpPtsPre)~=preNum && length(tmpPtsPost)~=postNum
                    continue;
                end
                [velocity, ~, totalTps] = velocityCal(tmpPtsPre, tmpPtsPost, Vec,...
                    curTps, j, Vec, curTps, j+1);
                %[~,ObzvVar] = obzVar(curTrack(j), tmpPtsPre, curTrack(j+1),...
                %    tmpPtsPost, totalTps, movieInfo);
                
                vCell{preNum+1, postNum+1}(cnt,1:3) = (Vec(j+1,:)-Vec(j,:))./(curTps(j+1)-curTps(j));
                vCell{preNum+1, postNum+1}(cnt,4:6) = velocity;
                vCell{preNum+1, postNum+1}(cnt,7:9) = [preNum, postNum, curTps(j+1)-curTps(j)];
                %vCell{preNum+1, postNum+1}(cnt,10:12) = ObzvVar;
                %vCell{preNum+1, postNum+1}(cnt,13:15) = curVar;
                cnt = cnt+1;
            end
        end
        vCell{preNum+1, postNum+1} = vCell{preNum+1, postNum+1}(1:cnt-1,:);
        
        v2v1Var(v2v1Cnt,1:4) = [preNum, postNum, 1, cnt-1];
        
        vv = vCell{preNum+1, postNum+1};
        v2v1Var(v2v1Cnt,5:7) = var(vv(:,1:3)-vv(:,4:6));
        totalTps = preNum+postNum;
        varRatio = (1+min(1,preNum)/totalTps)^2 + (min(1,preNum)/totalTps)^2+...
            (1+min(1,postNum)/totalTps)^2 + (min(1,postNum)/totalTps)^2;
        v2v1Var(v2v1Cnt,8) = varRatio;
        v2v1Var(v2v1Cnt,9:11) = v2v1Var(v2v1Cnt,5:7)./varRatio;
        v2v1Cnt = v2v1Cnt+1;
    end
end

idx = find(~isnan(v2v1Var(:,1)));
v0Var = v2v1Var(idx(2),5:7) - v2v1Var(idx(1),9:11)*v2v1Var(idx(2),8);
% v2v1Var= zeros(maxEdge^2,6);% preNum, postNum and jumpNum
% 
% for i=1:size(vvUnique,1)
%     zz = vv(:,7)==vvUnique(i,1) & vv(:,8)==vvUnique(i,2) & vv(:,9)==vvUnique(i,3);
%     v2v1Var(i,5:7) = var(vv(zz,1:3)-vv(zz,4:6));
%     v2v1Var(i,4) = sum(zz);
%     v2v1Var(i,8:10) = v2v1Var(i,5:7)-mean(vv(zz,10:12));
% end
[h,w]=size(vCell);
s = zeros(h*w,8);% prenum, postnum, k, b(y=kx*b)
ccMat = zeros(h,w);
for i=1:h
    for j=1:w
        if i==1 && j==1
            continue;
        end
        curCell = vCell{i,j};
        s((i-1)*w+j,1:2) = [i-1,j-1];
        for cc=1:3
            y = (curCell(:,cc)-curCell(:,cc+3)).^2;
            x = abs(curCell(:,cc+3));
            tmpKB = [x  ones(size(x))]\y;
            if cc==1
                ccRes = (corrcoef(x,y));
                ccMat(i,j) = ccRes(2);
            end
            %figure;scatter(x,y);hold on; 
            %plot(min(x):0.01:max(x), [min(x):0.01:max(x)]*tmpKB(1)+tmpKB(2));
            s((i-1)*w+j,1+2*cc:2+2*cc) = [tmpKB(1),tmpKB(2)];
        end
    end
end
end