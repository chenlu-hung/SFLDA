% manifold kernel smoothing, predict functional values/FPC scores based on
% functional manifold components
%
% y - functional manifold components
% S - a structure, see maniMDS
% kernel - smoothing kernel [default epan]
% h - smoothing bandwidth [default by Cross-Validation]
% out - leave-one-out indicator
% d - manifold dimensions used [default size(y,2)]
% score - indicator, use 1 if predict FPC scores
% x - predicted functional values/FPC scores

function x = maniKS(y,S,kernel,h,out,d,score)

if nargin<7|isempty(score) score=0; end
if nargin<6|isempty(d) d=size(y,2); end
if nargin<5|isempty(out) out=0; end
if nargin<4 h=[]; end
if nargin<3|isempty(kernel) kernel='epan'; end

d_mani = L2_distance(y(:,1:d)',S.Y(:,1:d)',0);
if isempty(h) h = bwCV(S,kernel); end
if isa(h,'char')
    K = min(str2num(h),S.N-out);
    if out == 1
        outnum = sum(d_mani'==0);
    else
        outnum = zeros(1,size(y,1));
    end
    [Sort_dis,Ord] = sort(d_mani');
    if K==1
        x = S.X_reg(Ord(1,:)+outnum,:);
        return;
    end
    H = Sort_dis(K+outnum,:);
else
    H = repmat(h,[1,size(y,1)]);
end

w = kernelval(repmat(1./H',[1,S.N]).*d_mani,kernel);
if out==1 w(find(w==kernelval(0,kernel)))=0; end
if score==1
    X = S.xi;
else
    X = S.X_reg;
end
if all(sum(w,2))
    x = repmat(1./sum(w,2),[1,S.N]).*w*X;
else
    idx = find(sum(w,2)>0);
    x = repmat(NaN,[size(y,1),S.M]);
    x(idx,:) = repmat(1./sum(w(idx,:),2),[1,S.N]).*w(idx,:)*X;
end

end