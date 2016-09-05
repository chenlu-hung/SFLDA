function fit = SFLDA(gnd,fea,tt,options)
% -------------------------------------------------------------------------
% Main program for sensible LDA
% Inputs
% gnd:           n*1 vector; groud truth labels
% fea:           1*n cell array; training trajectories
% tt:            1*n cell array; observatory times for training trajectories
% options:       structure array; see the examples for more details
% Output
% fit:           structure array; SLDA results
% -------------------------------------------------------------------------
    C = unique(gnd);
    nc = length(C);
    FVE_Threshold = options.FVE_threshold;
    q = options.q;

    [MU,FPCAObj] = doEstimation(gnd,fea,tt,options);
    fit_n = NSLDA(gnd,MU,FPCAObj);
    
    % q-fold cross validation
    if size(fit_n.eigvec,2) >= (nc-1)
        errors = qCV(q,gnd,fea,tt,options)
        if errors(1) < errors(2)
            fit = fit_n;
            disp('null LDA')
        else
            fit = DLDA(gnd,MU,FPCAObj);
            disp('direct LDA')
        end
    else
        % see (2.5)
        disp('combined')
        Gamma = zeros(size(MU));
        out21 = getVal(FPCAObj,'out21');
        dt = out21(2) - out21(1);
        noeig = getVal(FPCAObj,'no_opt');
        phi = getVal(FPCAObj,'eigen');
        for k = 1:nc
            Gamma(k,:) = MU(k,:) - MU(k,:)*fit_n.eigvec*fit_n.eigvec';
        end
        fit_d = DLDA(gnd,Gamma,FPCAObj);
        fit.eigvec = [fit_n.eigvec,fit_d.eigvec];
    end
    
    fit.labels = C;
    fit.mu = MU;
    fit.FPCAobj = FPCAObj;
    fit.fea = fea';
    fit.tt = tt';
    fit.gnd = gnd;
end
