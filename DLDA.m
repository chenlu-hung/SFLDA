function fit = DLDA(gnd,MU,FPCAobj)
% -------------------------------------------------------------------------
% Subroutine to perform LDA in the range space of within covariance operator
% Inputs
% gnd:           n*1 vector; groud truth labels
% MU:            c*p matrix; estimated residual functions
% FPCAobj:       PACE output for within covaraince function
% Output
% fit:           structure; results from (3.4)
% -------------------------------------------------------------------------
    labels = unique(gnd);
    [c,p] = size(MU);
    out21 = getVal(FPCAobj,'out21');
    dt = out21(2) - out21(1);
    Within = getVal(FPCAobj,'xcov');
    [UB,DB,VB] = svd(bsxfun(@minus,MU,mean(MU)),'econ');
    DB = diag(DB); DB = dt*(DB.^2); 
    idx = DB>1e-16; DB = DB(idx); VB = VB(:,idx);
    VB = VB / sqrt(dt);
    Omega_B = diag(DB);
    Omega_W = zeros(length(DB));
    for i = 1:length(DB)
        for j = 1:length(DB)
            Omega_W(i,j) = trapz(out21,VB(:,i).*(Within*VB(:,j)));
        end
    end
    [V,D] = eig(Omega_B,Omega_W);
    eigvec = VB*V;
    
    fit.eigvec = eigvec;
    fit.mu = MU;
    fit.labels = labels;
    fit.FPCAobj = FPCAobj;
end