clear all;

% addpath(genpath(path/to/PACE))

n = 100; p = 200;
[fea_Train,gnd_Train,tt_Train] = genSparse(n,p,3);
[fea_Test,gnd_Test,tt_Test] = genSparse(n,p,3);

options = setOptions('selection_k','FVE','FVE_threshold',0.95,'ngrid',p,'regular',0,'verbose','on','bwmu_gcv',0,'bwxcov_gcv',0);
options.q = 5;
fit = SFLDA(gnd_Train,fea_Train',tt_Train',options);
gnd_Predict = predict_SFLDA(fit,fea_Test,tt_Test);
sum(gnd_Predict~=gnd_Test)/length(gnd_Test)