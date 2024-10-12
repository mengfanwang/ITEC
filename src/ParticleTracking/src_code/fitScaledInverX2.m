function [df, tau2] = fitScaledInverX2(devXYZ)
% fit the scaled inverse chi-square distribution parameters
numDim = size(devXYZ,2);
numSample = size(devXYZ,1);
df = nan(numDim,1);
tau2 = nan(numDim,1);
for i = 1:numDim
    tau2(i) = numSample/sum(1./devXYZ(:,i));
    df(i) = mean(devXYZ(:,i))/(mean(devXYZ(:,i))-tau2(i));
    dfPre = inf;
    rhand = sum(log(devXYZ(:,i))) - numSample*log(tau2(i));
    while abs(df(i)-dfPre)>1e-6
        disp(abs(df(i)-dfPre));
        dfPre = df(i);
        df(i) = df(i)-(log(df(i))-psi(df(i))-rhand)/(1/df(i)-psi(1,df(i)));
    end
end
df = df*2;
end