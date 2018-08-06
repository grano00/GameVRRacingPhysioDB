function [ret] = myGaussianProcessRegression(XT,yT,xt,yt)
ops = struct('UseParallel',true);
tic
md = fitrgp(XT,yT,'KernelFunction','ardrationalquadratic');
ttime = toc;
%'HyperparameterOptimizationOptions',ops);
tic
pr = predict(md,xt);
ptime = toc;
err = sqrt(immse(yt,pr));
ret = {err,pr,yT,XT,yt,xt,md,ttime,ptime};
end