function [fea,gnd,tt] = genDense(n,p)
tgrid = linspace(0,1,p);
mu = zeros(3,p);


neig = 2;
phi = zeros(neig,p);

for i=1:neig
    phi(i,:) = sin(2*pi*tgrid*i);
end

fea = cell(3*n,1);
gnd = ones(3*n,1);
tt = cell(3*n,1);

% case c
 mu(1,:) = 0.2*cos(2*pi*tgrid);
 mu(2,:) = 0.2*cos(4*pi*tgrid);

% case b
% mu(1,:) = sin(2*pi*tgrid);
% mu(2,:) = sin(2*pi*tgrid)+0.25*cos(2*pi*tgrid);


% case a
% mu(1,:) = sin(2*pi*tgrid);
% mu(2,:) = sin(4*pi*tgrid);


for i = 1:n
    fea{i} = mu(1,:) +  normrnd(0,1/11,1,p);
    fea{n+i} = mu(2,:) + normrnd(0,1/11,1,p);
    fea{2*n+i} = mu(3,:) +  normrnd(0,1/11,1,p);
    
    for j=1:neig
        fea{i} = fea{i} + normrnd(0,1)/j*phi(j,:);
        fea{n+i} = fea{n+i} + normrnd(0,1)/j*phi(j,:);
        fea{2*n+i} = fea{2*n+i} + normrnd(0,1)/j*phi(j,:);
    end
    
    tt{i} = tgrid;
    tt{n+i} = tgrid;
    tt{2*n+i} = tgrid;
    gnd(i) = 1;
    gnd(n+i) = 2;
    gnd(2*n+i) = 3;
end