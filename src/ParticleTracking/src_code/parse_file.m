function [arcs, sst_init, distances, upt_n_nodes] = parse_file(arc_name, sst_name, dist_name, upt_sz_name)

% read arcs from original data
fileID = fopen(arc_name,'r');
for j=1:9
    tline = fgets(fileID);
end
arcs = fscanf(fileID, 'a %f %f %f\n');
arcs = reshape(arcs, 3, [])';
fclose(fileID);

% read initial shortest path tree
fileID = fopen(sst_name,'r');
sst = fscanf(fileID, '%f %f\n');
sst_init = sst(2:2:end);
fclose(fileID);

% read updated tree size in each iteration
fileID = fopen(upt_sz_name,'r');
upt_n_nodes = fscanf(fileID, '%f\n');
fclose(fileID);

% read all distances in sst
distances = cell(length(sst_init),1);
fileID = fopen(dist_name,'r');
distMat = fscanf(fileID, '%f, %f, %f\n');
fclose(fileID);

distMat = reshape(distMat, 3, [])';
cnt  = 0;
for i=1:2:max(distMat(:,2))
    val_id = distMat(:,2)==i;
    if isempty(find(val_id == 1,1))
        continue;
    end
    cnt = cnt + 1;
    distances{cnt} = distMat(val_id,3);
    distMat(val_id,:) = [];
end


st_len = cellfun(@length, distances);
[~, od] = sort(st_len, 'descend');
distances = distances(od);