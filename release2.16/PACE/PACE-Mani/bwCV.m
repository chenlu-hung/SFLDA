% 10 fold cross validation function to choose bandwidth h
%
% S - a structure, see maniMDS
% kernel - smoothing kernel
% bw - selected bandwidth

function bw = bwCV(S,kernel)

[h,Hcv] = h10cv(S,kernel);
[K,Kcv] = hknn10cv(S,kernel);
if min(Kcv{2})>min(Hcv{2})
    bw = h;
else
    bw = K;
end

end          
        