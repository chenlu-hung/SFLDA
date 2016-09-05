function fea_Reconstruct = Reconstruct(fit,fea_Test,tt_Test)
    NTest = length(fea_Test);
    fea_Reconstruct = cell(NTest,1);
    nc = length(fit.labels);
    FPCAobj = fit.FPCAobj;
    out21 = getVal(FPCAobj,'out21');
    xbar = fit.mu;
    sigma = getVal(FPCAobj,'sigma');
    lambda = getVal(FPCAobj,'lambda');
    phi = getVal(FPCAobj,'eigen');
    mu = mean(xbar,1);
    
    for i = 1:NTest
        logpseudolik = zeros(nc,1);
        tmp = cell(nc,1);
        fea_Reconstruct{i} = zeros(1,length(out21));
        tt = tt_Test{i};
        for j = 1:nc
            res{1} = fea_Test{i}-interp1(out21,xbar(j,:),tt,'spline');
            [junk,xi_new] = FPCApred(FPCAobj,res,tt);
            tmp{j} = xbar(j,:) + xi_new*phi';
            VW = interp1(out21,phi,tt,'spline');
            Within = VW*diag(lambda)*VW' + sigma*eye(length(tt));
            logpseudolik(j) = -res{1}*pinv(Within)*res{1}';
        end
        logpseudolik = logpseudolik - max(logpseudolik);
        pseudolik = exp(logpseudolik);
        for j = 1:nc
            fea_Reconstruct{i} = fea_Reconstruct{i} + pseudolik(j)*tmp{j}/sum(pseudolik);
        end
    end
end

