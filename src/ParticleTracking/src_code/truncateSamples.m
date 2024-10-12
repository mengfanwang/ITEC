function finalStd = truncateSamples(devXYZsingle, priorPara, g)
% the trajectory contains outlier edges with large value and lead to larger
% variance. We need to remove them before calculating variance
rmRatio = g.truncatedGaussian;
if isempty(priorPara) % if no gammaPrior input
    priorFlag = false; % for way 2 only
else
    priorFlag = true;
end
%% way 1 remove samples with value > rmRatio*std
if rmRatio > 1
    correctFactor = 1-rmRatio*normpdf(rmRatio)*2/(normcdf(rmRatio)-normcdf(-rmRatio));
    correctFactor = sqrt(correctFactor);
    finalStd = sqrt(mean(devXYZsingle.^2));
    tmpStd = [-1 -1 -1];
    devAbs = abs(devXYZsingle);
    while sum(abs(tmpStd-finalStd))>0.1
        tmpStd = finalStd;
        tmpThres = rmRatio*tmpStd;
        if sum(devAbs(:,1)>tmpThres(1))>0
            finalStd(1) = sqrt(mean(devAbs(devAbs(:,1)<=tmpThres(1),1).^2));%std(devAbs(devAbs(:,1)<=tmpThres(1),1));
        else
            finalStd(1) = sqrt(mean(devAbs(:,1).^2));
        end
        if sum(devAbs(:,2)>tmpThres(2))>0
            finalStd(2) = sqrt(mean(devAbs(devAbs(:,2)<=tmpThres(2),2).^2));%std(devAbs(devAbs(:,1)<=tmpThres(1),1));
        else
            finalStd(2) = sqrt(mean(devAbs(:,2).^2));
        end
        
        if sum(devAbs(:,3)>tmpThres(3))>0
            finalStd(3) = sqrt(mean(devAbs(devAbs(:,3)<=tmpThres(3),3).^2));%std(devAbs(devAbs(:,1)<=tmpThres(1),1));
        else
            finalStd(3) = sqrt(mean(devAbs(:,3).^2));
        end
        finalStd = finalStd./correctFactor;
    end
else
    %% way 2 remove rmRatio% samples
    % problem: for some trajectories, we do not need to remove particles
    numSample = size(devXYZsingle,1);
    if rmRatio*numSample<1
        numRm = 0;
        correctFactor = 1;
    else
        numRm = round(rmRatio*numSample);
        realRmRatio = (numRm/numSample+numRm/(numSample-1))/2;
        zValRt = abs(norminv(realRmRatio/2));
        zValLf = -zValRt;
        alpha = zValLf;
        beta = zValRt;
        correctFactor = 1 + (alpha*normpdf(alpha)-beta*normpdf(beta))/(normcdf(beta)-normcdf(alpha))...
            -((normpdf(alpha)-normpdf(beta))/(normcdf(beta)-normcdf(alpha)))^2;
    end
    finalStd = std(devXYZsingle);
    if priorFlag && size(devXYZsingle,1)>10 && finalStd(1)>3  
        %keyboard;
    end
    if priorFlag
        for i=1:3
            [~, stOd] = sort(abs(devXYZsingle(:,i)),'descend');
            valDev = devXYZsingle(stOd(numRm+1:end),i);
            tmpStd = stdFromPrior(valDev, priorPara(i,:),g);
            finalStd(i) = tmpStd/sqrt(correctFactor);
        end
    else
        if rmRatio*numSample<1
            return;
        end
        for i=1:3
            [~, stOd] = sort(abs(devXYZsingle(:,i)),'descend');
            valDev = devXYZsingle(stOd(numRm+1:end),i);
            tmpStd = sqrt((sum(valDev.^2)/(length(valDev)-1)));%std(valDev);%
            %alpha = -devXYZsingle(stOd(numRm+1),i)/tmpStd;%min(valDev)/tmpStd;%
            %beta = devXYZsingle(stOd(numRm+1),i)/tmpStd;%max(valDev)/tmpStd;%
            finalStd(i) = tmpStd/sqrt(correctFactor);
        end
    end
end
%% way 3: combine edge-cost change, standard deviation and percentage
% problem: for some trajectories, we do not need to remove particles