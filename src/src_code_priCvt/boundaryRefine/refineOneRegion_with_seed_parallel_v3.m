function [newLabel, comMaps] = refineOneRegion_with_seed_parallel_v3(seed_id, comMaps,option)
%% 20230616, Wei Zheng
%% gap detection based on min cut
comMaps.score3dMap(comMaps.score3dMap<0) = 0;
min_pv = min(comMaps.score3dMap(comMaps.score3dMap>0));
comMaps.score3dMap(comMaps.score3dMap<min_pv)=min_pv ;
max_pv = max(comMaps.score3dMap(:));
comMaps.score3dMap = scale_image(comMaps.score3dMap, 1e-3,1, min_pv, max_pv);

strel_rad = option.shift;
strel_ker = getKernel(strel_rad);
fMap = comMaps.regComp;
fMap = imdilate(fMap,strel_ker);
fMap((comMaps.idComp ~= seed_id) & (comMaps.idComp > 0)) = 0;
com = bwconncomp(comMaps.newIdComp);
newLabel = zeros(size(fMap));
for ii = 1:com.NumObjects
    sMap = false(size(fMap));
    sMap(com.PixelIdxList{ii}) = true;
    tMap = true(size(fMap));
    tMap(fMap) = false;
    tMap(comMaps.newIdComp) = true;
    tMap(com.PixelIdxList{ii}) = false;
    scoreMap = comMaps.score3dMap; % before fusion
%     scoreMap = comMaps.score3dMap; % after fusion
    scoreMap(scoreMap<0) = 0;
    [dat_in, src, sink] = graphCut_negHandle_mat(scoreMap, fMap, sMap, ...
        tMap, 10, [1 option.costDesign], true);       % before fusion
%     [dat_in, src, sink] = graphCut_negHandle_mat(scoreMap, fMap, sMap, ...
%         tMap, 26, [1 1], true);       % after fusion
    G = digraph(dat_in(:,1),dat_in(:,2),dat_in(:,3));
    if ~isempty(find(isnan(dat_in(:)), 1)) || isnan(sink)
        continue;
    end
    [~,~,cs,~] = maxflow(G, src, sink); % cs: fg, ct: bg
    cs = cs(cs<=numel(newLabel));
    %cs = cs(fMap(cs));
    cur_l = newLabel(cs);
    if length(unique(cur_l))>2
        keyboard;
    end
    newLabel(cs) = ii;
end


end

function strel_ker = getKernel(strel_rad)
    strel_ker = ones(strel_rad*2+1);
    ratio=size(strel_ker,1)/size(strel_ker,3);
    [xx,yy,zz] = ind2sub_direct(size(strel_ker), find(strel_ker));
    dist = sqrt( (xx - strel_rad(1)-1).^2 + (yy - strel_rad(2)-1).^2 + ( (zz -strel_rad(3)-1)*ratio ).^2 );
    strel_ker(dist>=strel_rad(1)) = 0;
    strel_ker = strel(strel_ker);
end