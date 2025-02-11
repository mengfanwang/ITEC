function movieInfo = mccTracker(dat_in, movieInfo, g, linkjumpoverCells, ...
    upt_cost_flag)
% min-cost circulation for solving min-cost flow problem
if nargin == 3
    linkjumpoverCells = false;
    upt_cost_flag = true; % default, after updating jumping ratio, we will update the edge cost
end
if true
    %% use CINDA framework
    if numel(dat_in) == 2
        % dat_in{1} = detection_arcs;
        % dat_in{2} = transition_arcs;
        [trajectories, costs] = mcc4mot(dat_in{1},dat_in{2});
        t_l = cellfun(@length, trajectories);
        trajectories(t_l<2) = [];
        costs(t_l<2) = [];
        ele1st = cellfun(@(x)x(1),trajectories);
        [~, od] = sort(ele1st,'ascend');
        trajectories = trajectories(od);
        costs = costs(od); % does not need indeed
        parents = cell(length(movieInfo.xCoord),1);
        kids = cell(length(movieInfo.xCoord),1);
        for i=1:numel(trajectories)
            % find the parent nodes and kid nodes for each cell
            for j=1:numel(trajectories{i})-1
                kids{trajectories{i}(j)} = cat(1, kids{trajectories{i}(j)}, trajectories{i}(j+1));
                parents{trajectories{i}(j+1)} = cat(1, parents{trajectories{i}(j+1)}, trajectories{i}(j));
            end
        end
        movieInfo.parents = parents;
        movieInfo.kids = kids;

    elseif numel(dat_in) == 3
        % dat_in{1} = detection_arcs;
        % dat_in{2} = transition_arcs;
        % dat_in{3}: links allows multi-cell to one
        [trajectories, costs, track_bf_merge, parents, kids] = ...
            mcc4mot_cell(dat_in{1},dat_in{2}, dat_in{3});
        t_l = cellfun(@length, trajectories);
        trajectories(t_l<2) = [];
        costs(t_l<2) = [];
        
        trajectories = cellfun(@(x) x(:), trajectories, ...
            'UniformOutput',false);
        movieInfo.track_bf_merge = track_bf_merge;
        movieInfo.parents = parents;
        movieInfo.kids = kids;
    end
    numTracks = numel(trajectories);
    % re-rank trajectories by their starting point
    
else
    %% shell based solver
    file_name = writeCirculationInputFile(dat_in, g);
    cd /home/congchao/Downloads/cs2-master;
    [status, cmdOut] = system(['./cs2 < ', file_name]);
    
    cd /media/congchao/32E20401E203C7D5/Users/Congchao/Desktop/Probjects_Google_Drive/ParticleTracking/;
    if status ~= 0
        error('cs2 function failed!\n');
    end
    
    [cost, dat1] = parsecs2OutputFile();
    
    %%%% backtrack tracks to get ids
    tmp   =  dat1(:, 1) == 1;
    start = dat1(tmp, 2);       %% starting nodes; is even
    
    tmp   =  mod(dat1(:, 1), 2) .* ~mod(dat1(:, 2), 2) .* (dat1(:, 2)-dat1(:, 1) ~= 1) ;
    links = dat1(tmp>0, 1:2);     %% links; is [even  odd]
    
    numTracks = length(start);
    trajectories = cell(numTracks,1);
    
    n_nodes = g.n_nodes;
    for i = 1:length(start)    %% for each track
        trajectories{i} = zeros(1e5,1);
        k = 0;
        this1 = start(i);
        while this1 ~= n_nodes
            k = k+1;
            trajectories{i}(k) = this1/2;
            this1 = links(links(:,1) == this1+1, 2);  %% should have only one value
            if (mod(this1, 2) + (length(this1) ~= 1) )>0         %% sanity check
                disp('error in the output of solver');
            end
        end
        trajectories{i} = trajectories{i}(1:k);
    end
end
if ~exist('particle2track','var')
    particle2track = nan(length(movieInfo.xCoord), 3); % 1: the track it belongs to; 2: its position 
    % in the track; 3: the original track before track merge it belongs to
    for i=1:numTracks
        particle2track(trajectories{i},1:2) = [i+zeros(length(trajectories{i}),1),...
            [1:length(trajectories{i})]'];
    end
end
if exist('track_bf_merge','var')
    for i=1:numel(track_bf_merge)
        particle2track(track_bf_merge{i},3) = i;
    end
end
movieInfo.tracks = trajectories;
movieInfo.pathCost = nan(numTracks,1);
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
if numel(dat_in) == 2 % linking allowing jump
    if g.loopCnt > g.maxIter + 2  % last two iter don't allow jump (Not used)
        ratio = [1 0 0];
    else
        ratio = zeros(1,g.k);
        for i=1:g.k
            ratio(i) = sum(gap==i)/length(gap);
        end
    end
    movieInfo.jumpRatio = ratio;
    if upt_cost_flag
        % update arc cost
        movieInfo = upt_cost_with_Dist(movieInfo);
        % the stable arcs
        movieInfo = stable_arc_cost_extract(movieInfo, g);
    end
    if linkjumpoverCells
        % for one to one linking, there may be un-necessary jump
        movieInfo = relinkJumpOveredCell(movieInfo, g);
    end
end
%g.jumpCost = movieInfo.ratio;
fprintf('finish one iteration of tracking!\n');

end