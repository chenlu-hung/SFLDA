function [xbar,FPCAObj] = doEstimation(gnd,fea,tt,options)
    C = unique(gnd);
    nc = length(C);

    if (options.regular==0) % Sparse data
        p = options.ngrid;
    else % dense data
        p = length(fea{1});
        options.ngrid = p;
    end

    xbar = zeros(nc,p);
    Res = cell(1,length(gnd));
    NTrain = zeros(1,nc);
    tt2 = cell(1,length(gnd));

    for i = 1:nc
        NTrain(i) = length(find(gnd==C(i)));
        y = fea(find(gnd==C(i)));
        t = tt(find(gnd==C(i)));
        %% calculate the mean for class i
        [xbar(i,:) out1 bw_mu] = smean(y,t,[],options.bwmu_gcv,options.ntest1,p,options.regular,[],0,options.verbose);
    end

    % Calculate Sigma_W
    k = 1;
    for i = 1:nc
        y = fea(find(gnd==C(i)));
        t = tt(find(gnd==C(i)));
        for j=1:NTrain(i)
            Res{k} = y{j} - interp1(out1,xbar(i,:),t{j},'spline','extrap');
            tt2{k} = t{j};
            k = k + 1;
        end
    end
    FPCAObj = myFPCA(Res,tt2,options);
end