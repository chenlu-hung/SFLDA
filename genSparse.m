function [fea,gnd,tt] = genSparse(n,p,Case)
tgrid = linspace(0,1,p);
fea = cell(3*n,1);
gnd = ones(3*n,1);
tt = cell(3*n,1);
mu = zeros(3,p);
switch Case
    case 1
        mu(1,:) = sin(2*pi*tgrid);
        mu(2,:) = sin(4*pi*tgrid);
    case 2
        mu(1,:) = sin(2*pi*tgrid);
        mu(2,:) = sin(2*pi*tgrid) + 0.25*cos(2*pi*tgrid);
    case 3
        mu(1,:) = 0.2*cos(2*pi*tgrid);
        mu(2,:) = 0.2*cos(4*pi*tgrid);
end

neig = 10;
for i = 1:n
    tmp = mu(1,:);
    for j = 1:neig
        tmp = tmp + normrnd(0,1)/j*sin(2*j*pi*tgrid);
    end
    tmp = tmp + normrnd(0,1,1,p)/11;
    ni = 1 + unidrnd(9);
    ts = sort(randperm(p,ni));
    fea{i} = tmp(ts);
    tt{i} = tgrid(ts);
    gnd(i) = 1;
    
    tmp = mu(2,:);
    for j = 1:neig
        tmp = tmp + normrnd(0,1)/j*sin(2*j*pi*tgrid);
    end
    tmp = tmp + normrnd(0,1,1,p)/11;
    ni = 1 + unidrnd(9);
    ts = sort(randperm(p,ni));
    fea{n+i} = tmp(ts);
    tt{n+i} = tgrid(ts);
    gnd(n+i) = 2;
    
    tmp = mu(3,:);
    for j = 1:neig
        tmp = tmp + normrnd(0,1)/j*sin(2*j*pi*tgrid);
    end
    tmp = tmp + normrnd(0,1,1,p)/11;
    ni = 1 + unidrnd(9);
    ts = sort(randperm(p,ni));
    fea{2*n+i} = tmp(ts);
    tt{2*n+i} = tgrid(ts);
    gnd(2*n+i) = 3;
end