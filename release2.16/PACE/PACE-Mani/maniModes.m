% function estimating manifold modes
%
% alpha - location coefficient (standardized); if 0, returns only
% the manifold mean
% S - a structure, see maniMDS
% fig - indicator, whether show plots
% kernel - smoothing kernel 
% h - smoothing bandwidth

function [modes,modes0,t] = maniModes(alpha,S,fig,kernel,h)

if nargin<5 h=S.h; end
if nargin<4 kernel=S.kernel; end
if nargin<3 fig=0; end

t = S.t;
if length(alpha)==1&alpha==0
    modes = maniKS(zeros(1,S.d),S,kernel,h,0);
    modes0 = S.mu;
    if fig==1
        figure
        hold on
        plot(S.t,modes0,'-k','LineWidth',2)
        plot(S.t,modes,'-b','LineWidth',2)
        legend('traditional mean','manifold mean','Location','Best')
        for k=1:length(S.Outliers)
            if ~S.Outliers(k)
                plot(S.T{k},S.X{k},'r-');
            end
        end
        plot(S.t,modes0,'-k','LineWidth',2)
        plot(S.t,modes,'-b','LineWidth',2)
        hold off
    end
else
    modes = {};
    modes0 = {};
    for j=1:S.d
        modes0{j} = repmat(S.mu,[length(alpha),1])+sqrt(S.lambda(j))*alpha'*S.phi(:,j)';
        for i=1:length(alpha)
            tmpt = zeros(1,S.d);
            tmpt(j) = alpha(i)*std(S.Y(:,j));
            modes{j}(i,:) = maniKS(tmpt,S,kernel,h,0);
            clear tmpt;
        end
    end
    if fig==1
        figure
        for k = 1:S.d
            subplot(S.d,2,2*k-1)
            plot(S.t,modes0{k})
            title(['variation due to FPC ' num2str(k)])
            subplot(S.d,2,2*k)
            plot(S.t,modes{k})
            title(['variation due to FMC ' num2str(k)])
        end
    end
end


end

