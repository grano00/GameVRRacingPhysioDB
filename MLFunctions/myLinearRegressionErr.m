function [ret] = myLinearRegressionErr(XT,yT,xt,yt)
ops = statset('UseParallel',true);
tic;
lr = fitrlinear(XT,yT,'Learner','leastsquares','Solver','sgd');
ttime = toc;
tic;
pr = predict(lr,xt);
ptime = toc;
err = sqrt(immse(yt,pr));
ret = {err,pr,yT,XT,yt,xt,lr,ttime,ptime};
end