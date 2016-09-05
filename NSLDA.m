function fit = NSLDA(gnd,MU,FPCAobj)
% -------------------------------------------------------------------------
% Subroutine to perform LDA in the null space of within covariance operator
% Inputs
% gnd:           n*1 vector; groud truth labels
% MU:            c*p matrix; estimated mean functions
% FPCAobj:       PACE output for within covaraince function
% Output
% fit:           structure; results from (3.3)
% -------------------------------------------------------------------------
    labels = unique(gnd);
    [c,p] = size(MU);
    % see (2.3)
    options = getVal(FPCAobj,'ops');
    FVE_Threshold = options.FVE_threshold;
    out21 = getVal(FPCAobj,'out21');
    dt = out21(2) - out21(1);
    noeig = getVal(FPCAobj,'no_opt');
    lambda = getVal(FPCAobj,'lambda');
    phi = getVal(FPCAobj,'eigen');
    Gamma = MU;
    for k = 1:c
        for i = 1:noeig
            Gamma(k,:) = Gamma(k,:) - trapz(out21,MU(k,:).*phi(:,i)')*phi(:,i)';
        end    
    end
    % see (2.5)
    [U,DG,VG] = svd(Gamma,'econ');
    DG = diag(DG); 
    DG = dt*(DG.^2); idx = DG>1e-8; DG = DG(idx);
    VG = VG(:,idx)/sqrt(dt);
    FVE = cumsum(DG)./sum(DG);
    q = min([find(FVE>FVE_Threshold,1,'first'),c-1]);
    VG = VG(:,1:q); DG = DG(1:q);
    fit.eigvec = VG;
    fit.mu = MU;
    fit.labels = labels;
    fit.r = Gamma;
    fit.FPCAobj = FPCAobj;
end
