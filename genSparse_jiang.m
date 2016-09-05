function [fea,gnd,tt] = genSparse(n,p)
tgrid = linspace(0,1,p);
mu = zeros(3,p);

phi = zeros(10,p);
for i=1:10
    phi(i,:) = sin(2*pi*tgrid*i);
end

fea = cell(3*n,1);
gnd = ones(3*n,1);
tt = cell(3*n,1);

mu(1,:) = 0.2*cos(2*pi*tgrid);
mu(2,:) = 0.2*cos(4*pi*tgrid);

% mu(1,:) = sin(2*pi*tgrid);
% mu(2,:) = sin(2*pi*tgrid)+0.25*cos(2*pi*tgrid);

for i = 1:n
    tmp1 = mu(1,:) +  normrnd(0,1/11,1,p);
    tmp2 = mu(2,:) +  normrnd(0,1/11,1,p);
    tmp3 = mu(3,:) +  normrnd(0,1/11,1,p);
    
    for j=1:10
        tmp1 = tmp1 + normrnd(0,1/j,1,1)*phi(j,:);
        tmp2 = tmp2 + normrnd(0,1/j,1,1)*phi(j,:);
        tmp3 = tmp3 + normrnd(0,1/j,1,1)*phi(j,:);
    end
    
    ni1 = 2 + unidrnd(10); ts1 = sort(randperm(p,ni1));
    ni2 = 2 + unidrnd(10); ts2 = sort(randperm(p,ni2));
    ni3 = 2 + unidrnd(10); ts3 = sort(randperm(p,ni3));    
    
    fea{i} = tmp1(ts1);
    fea{n+i} = tmp2(ts2);
    fea{2*n+i} = tmp3(ts3);
    
    tt{i} = tgrid(ts1);
    tt{n+i} = tgrid(ts2);
    tt{2*n+i} = tgrid(ts3);
    gnd(i) = 1;
    gnd(n+i) = 2;
    gnd(2*n+i) = 3;
end