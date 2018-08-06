function [ret] = mySVMGaussianRegressionErr(XT,yT,xt,yt)
ops = statset('UseParallel',true);
tic;
svmg = fitrsvm(XT,yT,'KernelFunction','gaussian','KernelScale','auto')
ttime = toc;
tic;
pr = predict(svmg,xt);
ptime = toc;
err = sqrt(immse(yt,pr));
ret = {err,pr,yT,XT,yt,xt,svmg,ttime,ptime};
end