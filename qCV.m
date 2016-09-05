function errors = qCV(q,gnd,fea,tt,options)
% -------------------------------------------------------------------------
% q-fold cross validation tp select whether NSLDA or DLD should be used
% Inputs
% q:             scalar; q-fold
% gnd:           n*1 vector; groud truth labels
% fea:           1*n cell array; training trajectories
% tt:            1*n cell array; observatory times for training trajectories
% options:       structure array; see the examples for more details
% Outputs
% errors:        2*1 vector; error rates for NSLDA errors(1) and DLDA errors(2)
% -------------------------------------------------------------------------
    errors = zeros(1,2);
    part = cvpartition(gnd,'kfold',q);
    for i = 1:q
        gnd_Train = gnd(training(part,i));
        fea_Train = fea(training(part,i));
        tt_Train = tt(training(part,i));
        gnd_Test = gnd(test(part,i));
        fea_Test = fea(test(part,i));
        tt_Test = tt(test(part,i));
        [MU,FPCAObj] = doEstimation(gnd_Train,fea_Train,tt_Train,options);
        fit_n = NSLDA(gnd_Train,MU,FPCAObj);
        fit_n.fea = fea'; fit_n.tt = tt'; fit_n.gnd = gnd;
        errors(1) = errors(1) + sum(gnd_Test~=predict_SFLDA(fit_n,fea_Test',tt_Test))/length(gnd_Test);
        fit_d = DLDA(gnd_Train,MU,FPCAObj);
        fit_d.fea = fea'; fit_d.tt = tt'; fit_d.gnd = gnd;
        errors(2) = errors(2) + sum(gnd_Test~=predict_SFLDA(fit_d,fea_Test',tt_Test))/length(gnd_Test);
    end
    errors = errors/q;
end