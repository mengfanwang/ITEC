function [zScore,t0] = getOrderStatEasy( d, smask, bmask, tmask, qVar, qNtry, pMu, pSigma, pTest, pMinIntensityDiff )
%GETORDERSTAT Compute order statistics
% d: transformed data. smask: synapse mask.
% bmask: background mask, MUST be connected. imgVar: image variance

% find balanced groups -----
g1 = d(smask>0);
nPix = sum(smask(:));
tmask1 = tmask - smask;
tmask1(tmask1<0) = 0;
%nmask = (bmask - smask).*(1-tmask1);
% nNeibTot = sum(nmask(:));
% if nNeibTot < 4
%     zScore = 0;
%     return
% end
%nPix = 25;
if nPix < 4
    nPix = 4;  % !!
end
if qNtry>5
    nPix = qNtry;
end;
%% variance estimation -----
% X0 = smask;
% done = 0;
% ntry = 1;
% nNeibPrev = 0;
%tic;
% while ~done
%     X0 = imdilate(X0,strel('square',3));
%     X = X0.*(1-smask).*bmask.*(1-tmask1);
%     nNeib = sum(X(:));
% %     if nNeib>=nPix && nNeib>50
%     if nNeib>=nPix
%         g2 = d(X>0);
%         done = 1;
%     elseif (nNeibPrev==nNeib) || ntry==qNtry
%         g2 = d(X>0);
%         done = 1;
%     else
%         done = 0;
%     end
%     if done
%         g2 = d(X>0);
%     end
%     ntry = ntry + 1;
%     nNeibPrev = nNeib;
% end
% toc;
%%version 2 for variance estimation ----- This is much quicker

[Ny,Nx] = size(smask);
[Y0,X0] = find(smask>0);
Coord0 = [Y0,X0];
NeiCoord = [-1 -1 -1 0 -1 1 0 -1 0 1 1 -1 1 0 1 1];

done = 0;
ntry = 1;
nNeibPrev = 0;

while ~done
    % generate the dilated part of synapse
    Coord = repmat(Coord0,1,8);
    Coord = bsxfun(@plus,Coord,NeiCoord);
    Coord = [Coord0;Coord(:,1:2);Coord(:,3:4);Coord(:,5:6);Coord(:,7:8);...
        Coord(:,9:10);Coord(:,11:12);Coord(:,13:14);Coord(:,15:16)];
    OKCoord = Coord(:,1)>0 & Coord(:,1)<=Ny & Coord(:,2)<=Nx & Coord(:,2)>0;
    Coord0 = Coord(OKCoord,:);
    Coord0 = unique(Coord0,'rows');
    
    Coord_ind = Coord0(:,1)+(Coord0(:,2)-1)*Ny;
    
    X = (1-smask(Coord_ind)).*bmask(Coord_ind).*(1-tmask1(Coord_ind));
    nNeib = sum(X(:));
%     if nNeib>=nPix && nNeib>50
    if nNeib>=nPix
        g2 = d(Coord_ind(X>0));
        done = 1;
    elseif (nNeibPrev==nNeib) || ntry==qNtry
        g2 = d(Coord_ind(X>0));
        done = 1;
    else
        done = 0;
    end
    if done
        g2 = d(Coord_ind(X>0));
    end
    ntry = ntry + 1;
    nNeibPrev = nNeib;
end

%%
imgVar = qVar;
% imgVar = min(qVar,var(g2));
% imgVar = max(imgVar,0.005);
% if thr>0
%     g1=g1(g1>thr);
%     g2 = g2(g2<thr);
% end;
%% test statistics -----
if strcmp(pTest,'order')
    t0 = mean(g1) - mean(g2);
    M = length(g1);
    N = length(g2);
    if pMinIntensityDiff>0 && t0<pMinIntensityDiff
        zScore = 0;
        return
    end
    if M < 4
        M = 4;
    end
    
    if N<=3
        zScore = 0;
        return
    end
    
    if N < (M/10)
        zScore = 0;
        return
    end
    
    if M>100
        M = 100;
    end
    if N>100
        N = 100;
    end
    
    mu = pMu(M,N);
    sigma = pSigma(M,N);
    mu = mu*sqrt(imgVar);
    sigma = sigma*sqrt(imgVar);
    
    % [mu, sigma] = getOrderStatSingle( M, N, imgVar );
    zScore = (t0-mu)/sigma;
%     if zScore <=0
%         pvalue = inf;
%     else
%         pvalue = 1-normcdf(zScore,mu,sigma);
%     end;
    t0 = t0/sqrt(imgVar);
else
    if length(g1)>100
        idx = randperm(length(g1),100);
        g1 = g1(idx);
    end
    if length(g2)>100
        idx = randperm(length(g2),100);
        g2 = g2(idx);
    end
    [~,~,~,ss] = ttest2(g1,g2);
    zScore = ss.tstat;
end



end







