clear all;


nc = 2; % number of classes
m = [200 250]; % number of samples in each classes

p = 50; % number of observation per subject for functional case

alldata = cell(nc,1);
alltime = cell(nc,1);
%% generating data
%{
t = linspace(0,1,p);
y1 = cell(1,m(1));
t1 = cell(1,m(1));
ni = unidrnd(10,[1 m(1)]);
for i = 1:m(1)
    tmp = 0.5*cos(2*pi*t) + normrnd(0,1)*sin(2*pi*t) + normrnd(0,1/2)*sin(4*pi*t) + normrnd(0,1,1,p);
    ts = sort(randperm(p,ni(i)));
    y1{i} = tmp(ts);
    t1{i} = t(ts);
end

y2 = cell(1,m(2));
t2 = cell(1,m(2));
ni = unidrnd(10,[1 m(2)]);
for i = 1:m(2)
    tmp = normrnd(0,1)*sin(2*pi*t) + normrnd(0,1/2)*sin(4*pi*t) + normrnd(0,1,1,p);
    ts = sort(randperm(p,ni(i)));
    y2{i} = tmp(ts);
    t2{i} = t(ts);
end
%}
[fea,gnd,tt] = Jiang14Sparse([200 250],50);
y1 = fea(1:m(1)); t1 = tt(1:m(1));
y2 = fea((m(1)+1):(m(1)+m(2))); t2 = tt((m(1)+1):(m(1)+m(2)));

xbar = zeros(nc,p);
k=1;
Res = cell(1,sum(m));
tt = cell(1,sum(m));


%% calculate the mean for class 1
    [xbar(1,:) out1] = smean(y1,t1,[],[],[],p,0,[],[],[]);

% get the residuals for Sigma_W    
    for j=1:m(1)
        Res{k} = y1{j} - interp1(out1,xbar(1,:),t1{j});
        tt{k} = t1{j};
        k=k+1;
    end

%% calculate the mean for class 2   
    [xbar(2,:) out1] = smean(y2,t2,[],[],[],p,0,[],[],[]);  
    
% get the residuals for Sigma_W    
    for j=1:m(2)
        Res{k} = y2{j} - interp1(out1,xbar(2,:),t2{j});
        tt{k} = t2{j};
        k=k+1;
    end

Sigma_B = zeros(p);

mu = mean(xbar,1); % calculate the overall mean

%% calculate Sigma_B
for i=1:nc
    tmp = xbar(i,:)-mu;
    Sigma_B = Sigma_B + m(i)*(tmp'*tmp);
end

%% calculate Sigma_W
% [Sigma_W out21 Sigma_Raw] = SigW(Res,tt,[],[],[],p,0,1,[],[],[]);
options = setOptions(); options.ngrid = 50; options.regular = 0; options.error = 1;
fit = myFPCA(Res,tt,options);
% Sigma_W = getVal(fit,'xcov');

%{
figure;
surf(Sigma_B)

figure;
surf(Sigma_W)
%}
