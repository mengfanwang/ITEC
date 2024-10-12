function movieInfo = mccTracker(dat_in, movieInfo, g)
% min-cost circulation for solving min-cost flow problem

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
particle2track = nan(g.particleNum, 2); % the track it belongs to and position

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

for i=1:numTracks
    particle2track(trajectories{i},:) = [i+zeros(length(trajectories{i}),1),...
        [1:length(trajectories{i})]'];
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
ratio = zeros(1,g.k);
for i=1:g.k
    ratio(i) = sum(gap==i)/length(gap);
end
movieInfo.jumpRatio = ratio;
%g.jumpCost = movieInfo.ratio;
fprintf('finish one iteration of tracking!\n');
end