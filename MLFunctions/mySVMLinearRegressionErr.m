function [ret] = mySVMLinearRegressionErr(XT,yT,xt,yt)
ops = statset('UseParallel',true);
tic;
svml = fitrsvm(XT,yT,'KernelFunction','linear','KernelScale','auto')
ttime = toc;
tic;
pr = predict(svml,xt);
ptime = toc;
err = sqrt(immse(yt,pr));
ret = {err,pr,yT,XT,yt,xt,svml,ttime,ptime};
end