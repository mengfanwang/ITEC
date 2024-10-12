function estimatedStd = stdFromPrior(devXYZ, priorFuncPara, g)
% using prior distribution to estimate the standard deviation of a
% trajectory
if strcmp(g.priorType, 'gamma')
    % we use gamma distribution now, mean is 0
    dim = size(devXYZ,2);
    tau = nan(1,dim);
    for i=1:dim % for x,y,z directions
        x = devXYZ(:,i);
        alpha = priorFuncPara(i,1) + length(x)/2;
        beta = priorFuncPara(i,2) + 0.5*sum(x.^2);
        tau(i) = (alpha-1)/beta; % tau is precision=1/sigma^2
    end
    estimatedStd = sqrt(1./tau);
elseif strcmp(g.priorType, 'sInvX2')
    % we use scale-inverse chi-square distribution now, mean is 0
    dim = size(devXYZ,2);
    estimatedStd = nan(1,dim);
    for i=1:dim % for x,y,z directions
        x = devXYZ(:,i);
        newdf = priorFuncPara(i,1) + length(x);
        newTau2 = (priorFuncPara(i,1)*priorFuncPara(i,2) + sum(x.^2))/newdf;
        estimatedStd(i) = sqrt(newdf*newTau2/(newdf+2));
    end
    
elseif strcmp(g.priorType, 'weiAvg')
    dim = size(devXYZ,2);
    estimatedStd = nan(1,dim);
    for i=1:dim % for x,y,z directions
        x = devXYZ(:,i);
        estVmean = sum(x.^2)/(length(x)-1);
        estVvar = estVmean^2*2*length(x)/((length(x)-1)^2);
        wEstVmean = (priorFuncPara(i,2))/(priorFuncPara(i,2) + estVvar);
        wPrior = 1-wEstVmean;
        estimatedStd(i) = sqrt(wPrior*(priorFuncPara(i,1))+wEstVmean*estVmean);
    end
elseif strcmp(g.priorType, 'Gauss')
    dim = size(devXYZ,2);
    estimatedStd = nan(1,dim);
    for i=1:dim % for x,y,z directions
        x = devXYZ(:,i);
        n = length(x);
        sigma = sqrt(sum(x.^2)/(n-1)); % for sigma, we 
        estStdMean = sqrt(2/n)*gamma(0.5*n)*sigma/gamma(0.5*(n-1));
        estStdVar = g.dependencyFactor*(1/n)*(n-1-2*gamma(0.5*n)^2/gamma(0.5*(n-1))^2)*sigma^2;
        wEstVmean = (priorFuncPara(i,2))/(priorFuncPara(i,2) + estStdVar);
        wPrior = 1-wEstVmean;
        estimatedStd(i) = wPrior*(priorFuncPara(i,1))+wEstVmean*estStdMean;
    end
elseif strcmp(g.priorType, 'GaussGauss')
% have not done
else
    error('Undefined prior!\n');
end
end