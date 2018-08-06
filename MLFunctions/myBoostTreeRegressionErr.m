function [ret] = myBoostTreeRegressionErr(XT,yT,xt,yt)
ops = statset('UseParallel',true);
t = templateTree('Surrogate','On');
tic;
bt = fitensemble(XT,yT,'LSBoost',100,t);
ttime = toc;
tic;
pr = predict(bt,xt);
ptime = toc;
err = sqrt(immse(yt,pr));
ret = {err,pr,yT,XT,yt,xt,bt,ttime,ptime};
end