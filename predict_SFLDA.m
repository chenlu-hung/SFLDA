function [gnd_Test] = predict_SFLDA(fit,fea_Test,tt_Test)
    fea_Train = fit.fea; tt_Train = fit.tt;
    if ( getVal(fit.FPCAobj,'regular')<2 )
        fea_Train = Reconstruct(fit,fea_Train,tt_Train);
        fea_Test = Reconstruct(fit,fea_Test,tt_Test);
    end
    fea_Train = cell2mat(fea_Train);
    fea_Test = cell2mat(fea_Test);
    out21 = getVal(fit.FPCAobj,'out21');
    dt = out21(2) - out21(1);
    [n,p] = size(fea_Test);
    gnd_Test = ones(n,1);
    K = size(fit.mu,1);
    noeig = size(fit.eigvec,2);

    newmu = zeros(K,noeig);
    for k = 1:K
        for j = 1:noeig
            newmu(k,j) = trapz(out21,fit.mu(k,:).*fit.eigvec(:,j)');
        end
    end
    newfea_Train = zeros(size(fea_Train,1),noeig);
    for i = 1:size(fea_Train,1)
        for j = 1:noeig
            newfea_Train(i,j) = trapz(out21,fea_Train(i,:).*fit.eigvec(:,j)');
        end
    end
    newfea_Test = zeros(n,noeig);
    for i = 1:n
        for j = 1:noeig
            newfea_Test(i,j) = trapz(out21,fea_Test(i,:).*fit.eigvec(:,j)');
        end
    end
    S = cov(newfea_Train);
    for i = 1:n
        dist = zeros(K,1);
        for k = 1:K
            dist(k) = (newfea_Test(i,:)-newmu(k,:))*pinv(S)*(newfea_Test(i,:)-newmu(k,:))';
        end
        [trash,idx] = min(dist);
        gnd_Test(i) = fit.labels(idx);
    end

end