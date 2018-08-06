function [ret] = myRandomForestRegressionErr(XT,yT,xt,yt)
ops = statset('UseParallel',true);
tic;
tb = TreeBagger(128,XT,yT,'Method','regression','Options',ops);
ttime = toc;
tic;
pr = predict(tb,xt);
ptime = toc;
err = sqrt(immse(yt,pr));
ret = {err,pr,yT,XT,yt,xt,tb,ttime,ptime};
end