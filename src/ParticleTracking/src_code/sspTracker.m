function [movieInfo, res_G] = sspTracker(orgG, movieInfo, g)
% successive shortest path for min-cost flow
vNum = g.particleNum;
res_G = orgG;
pSet = cell(vNum,1);
numTracks = 0;
tic;
for i=1:vNum
%     if i==940
%         keyboard;
%     end
    [P,d] = shortestpath(res_G,g.excess_node(1),g.excess_node(2)); % 1 is source, n_nodes is sink saved in g.excess_node
    if d>=-0
        break;
    end
%     if sum(P==14082*2+1)>0
%         %keyboard;
%         pp = P(mod(P,2)==0);
%         pp = pp/2;
%     end
    fprintf('The %d paths found with cost: %.2f\n', i, d);
    % change weight as -1*weight
%     for j=1:length(P)-1
%         tmpIdx = find(res_G.Edges.EndNodes(:,1) == P(j) & res_G.Edges.EndNodes(:,2) == P(j+1)); % should only have one
%         if length(tmpIdx) ~= 1
%             error('There should be only one edge be changed.');
%         end
%         res_G.Edges.Weight(tmpIdx) = -1*res_G.Edges.Weight(tmpIdx);
%     end
    edgeIdx = findedge(res_G, P(1:end-1),P(2:end));
    res_G.Edges.Weight(edgeIdx) = -1*res_G.Edges.Weight(edgeIdx);
    % inverse edge
    res_G = flipedge(res_G,P(1:end-1),P(2:end));
    numTracks = numTracks + 1;
    pSet{i} = P;
end
timeLapse = toc;
fprintf('finish successive shortest path alg with %d paths, %.3fs!\n',numTracks, timeLapse);

% recover the tracks
trajectories = cell(numTracks,1);    
stNodes = res_G.Edges.EndNodes(res_G.Edges.EndNodes(:,1)==g.n_nodes,2);
%SSP = cell(numTracks,1);
pathCost = zeros(numTracks,1);
particle2track = nan(vNum,2); % the track it belongs to and position
%stNodes = stNodes-1;
parfor i=1:numTracks
    head = stNodes(i);
    trajectories{i} = stNodes(i)-1;
    pathCost(i) = pathCost(i)+res_G.Edges.Weight(res_G.Edges.EndNodes(:,1) == head-1);
    head = res_G.Edges.EndNodes(res_G.Edges.EndNodes(:,1) == head-1, 2);
    if head == 1
        pathCost(i) = pathCost(i)-g.c_en;
    end
    while head > 1
        if length(head) ~= 1
            error('There should be only one edge.');
        end
        if res_G.Edges.EndNodes(res_G.Edges.EndNodes(:,1) == head-1, 2) ~= 1
            pathCost(i) = pathCost(i)+res_G.Edges.Weight(res_G.Edges.EndNodes(:,1) == head-1);
        end
        trajectories{i} = [trajectories{i}, head - 1];
        head = res_G.Edges.EndNodes(res_G.Edges.EndNodes(:,1) == head-1, 2);
    end

    trajectories{i}  = trajectories{i}(end:-1:1)/2;
%     if sum(trajectories{i}==276)
%         disp(i);
%     end
    pathCost(i) = -pathCost(i)+sum(movieInfo.Ci(trajectories{i}))+g.c_en+g.c_ex;
end

for i=1:numTracks
    particle2track(trajectories{i},:) = [i+zeros(length(trajectories{i}),1),...
        [1:length(trajectories{i})]'];
end
movieInfo.tracks = trajectories;
movieInfo.pathCost = pathCost;
movieInfo.particle2track = particle2track;

% update the jump cost
gap = cell(numel(movieInfo.tracks),1);
for i=1:numel(movieInfo.tracks)
    if length(movieInfo.tracks{i})<g.trackLength4var
        continue;
    end
    fms = movieInfo.frames(movieInfo.tracks{i});
    gap{i} = fms(2:end)-fms(1:end-1);
    gap{i} = gap{i}(:);
end
gap = cat(1,gap{:});
ratio = zeros(1,g.k);
for i=1:g.k
    ratio(i) = sum(gap==i)/length(gap);
end
movieInfo.jumpRatio = ratio;
%g.jumpCost = movieInfo.ratio;
fprintf('finish one iteration of tracking!\n');

end