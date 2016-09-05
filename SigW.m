
%================
%Input Arguments:
%================
%
%Input y:          1*n cell array, y{i} is the vector of measurements for the
%                  ith subject, i=1,...,n.                     
%     
%Input t:          1*n cell array, t{i} is the vector of time points for the
%                  ith subject for which corresponding measurements y{i} are
%                  available, i=1,...,n.
%
%Input bwxcov:     1*2 vector, bandwidths for covariance surface used for
%                  smoothing of cov(X(t),X(s))
%                  bwxcov(i): ith coordinate of bandwidth vector, i=1,2.
%                  bwxcov(1)==0 & bwxcov(2)==0: use cross-validation (CV)
%                  or generalized cross-validation (GCV) for automatic
%                  selection.                                  [Default]
%                  bwxcov(1)>0 & bwxcov(2)>0: user-specified bandwidths.
%
%Input bwxcov_gcv: If setting bwxcov = [0 0], automatic bandwidth selection
%                  by CV or GCV choices
%                  0: CV method (may be time-consuming, can be accelerated by 
%                     choosing small values for ntest1 and ngrid1).
%                  1: GCV method (faster)                      [Default]
%                  2: Geometric mean between the minimum bandwidth and GCV
%                     bandwidth
%
%Input ntest1:     integer(<=n), number of curves used for CV when choosing
%                  bandwidths for smoothing the covariance surface. The subjects 
%                  in the test set are randomly selected from n subjects. Small 
%                  ntest1 will accelerate CV at less accuracy. [Default is 30.]
%
%Input ngrid1:     integer, number of support points for the covariance surface 
%                  in the CV procedure (selecting bandwidths of covariance).
%                  Note that input ntest1 and ngrid1 provide options to save
%                  computing time when using CV or GCV.        [Default is 30.]

%
%Input regular:   0, sparse (or irregular) functional data.      
%                 1, regular data with missing values
%                 2, completely balanced (regular) data.
%                 [], automatically chosen based on the design in t.   [Default]
%
%Input error:     0, no additional measurement error assumed.
%                 1, additional measurement error is assumed.    [Default]
%
%
%Input kernel:    a character string to define the kernel to be used in the
%                 1-D or 2-D smoothing
%                 kernel = 'epan'  ==> Epanechnikov kernel [Default for dense designs with n_i >= 20]
%                          'rect'  ==> Rectangular kernel
%                          'gauss'  ==> Gaussian kernel    [Default for sparse designs, regular designs with
%                                                           missings, dense designs for n_i < 20]
%                 Note: The Gaussian kernel is overall best for sparse designs but is slower than the other kernels 
%                       and if computational speed is of concern then one may wish to use the Epanechnikov kernel 
%                       also for the case of sparse designs.
%
%Input numBins:   0: no binning
%                 a positive interger (>= 10): prebin the data with user-defined
%                 number of bins. When numBins < 10, no binning will be performed.
%		  []:  prebin the data with the following rule    [Default]
%
%                 i) When the input data is regular = 1 or 2
%                    m = max of n_i, where n_i is the number of repeated measurements
%                    for ith subject.
%                 ii) regular = 0
%                    m = median of n_i
%
%                 When m <= 20 subjects, no binning.
%                 When n <= 5000 subjects and m <= 400, no binning.
%                 When n <= 5000 subjects and m > 400, numBins = 400.
%                 When n > 5000 subjects, compute
%
%                 m* = max(20, (5000-n)*19/2250+400)
%
%                 if m > m*, numBins = m*
%                 if m <= m*, no binning
%
%                 This option accelerates analysis, especially for data with
%                 many time measurements.
%

%Input verbose:   a character string
%                 'on': display diagnostic messages       [Default]
%                 'off': suppress diagnostic messages
%
%

% =================
% Output Arguments:
% =================
% Output xcov: smoothed covariance matrix
% Output out21: evaluation grid for xcov
% Output rcov: raw covariance matrix


function [xcov,out21,rcov] = SigW(y,t,bwxcov,bwxcov_gcv,ntest1,ngrid,regular,error,kernel,numBins,verbose)

                xcov=[];

invalid = 0;
userrange = [];
[invalid_data]  = CheckData(y, t); %NaN values and number of subjects.
if (invalid_data == 1)
    return
end

%Set default values for the major input arguments
%when the following input arguments are set to 
%be "[]" by the user


if isempty(bwxcov)
    bwxcov = [0 0]; %bandwidth choices for covariance function is CV or GCV
end
if isempty(bwxcov_gcv)
    bwxcov_gcv = 1; %bandwidth choices for covariance function is GCV
end
if isempty(ntest1)
    ntest1 = 30;
end
if isempty(ngrid)
   ngrid = 30;
end


if isempty(error)
    error = 1;       %error assumption with measurement error
end


if isempty(verbose)
  verbose = 'on';
end

 ncohort=length(t);     % obtain the number of curves or subjects

 ni = zeros(1,ncohort);
 for i = 1:ncohort
   ni(i)= length(t{i});
 end
 if all(ni == 1)
     fprintf(1,'Error: FPCA is aborted because the data do not contain repeated measurements!');
     return;
 end

 if isempty(regular)||~any(regular==0:2)
        regular = isregular(t);
 else
        ireg = isregular(t);
        if ireg < regular
            switch ireg
                case 0
                    fprintf(1,'Warning: the design is sparse but has been specified as regular or regular with missing.  No computation is performed.  Please rerun with default design setting.\n');
                case 1
                    fprintf(1,'Warning: the design is regular with missing but has been specified as regular.  No computation is performed.  Please rerun with default design setting.\n');
            end
            return;
        end
 end


%Prebin the data if conditions are satisfied as specified in the help
%for numBins.
  
if isempty(numBins)
   [newy,newt] = binData(y,t,regular,verbose);
elseif numBins >= 10
   [newy,newt] = binData(y,t,regular,verbose,numBins);
elseif numBins == 0
   newy = [];   %no binning set by user
elseif numBins < 10 || numBins < 0
   newy = [];   %no binning due to number of bins is too small
   fprintf(1,'Warning: number of bins must be at least 10! No binning will be performed!\n');
end

if ~isempty(newy)
    y = newy;
    t = newt;
    if regular == 0
        regular = 1;
    end
end

if isempty(kernel)
    if regular == 2 && length(t{1}) >= 20
       kernel = 'epan'; %kernel: Epanechnikov
    else
       kernel = 'gauss';   %kernel: Gaussian
    end
else
    kernNames = {'rect','gauss','epan','gausvar','quar'};
    if isempty(strmatch(kernel, kernNames, 'exact'))
        fprintf(1,['Warning: kernel ' kernel ' is unrecognizable! Reset to default now.\n']);
        if regular == 2 && length(t{1}) >= 20
            kernel = 'epan';
        else
            kernel = 'gauss';
        end
    end
end



%pool all the subjects and their corresponding time points into 1 x N vectors
%tt: vector to hold time points
%yy: vector to hold the observed measurements
%ind: indices for each subject
tt = cell2mat(t);  % 1 x N vector to hold the observed time points from all subjects
yy = cell2mat(y);  % 1 x N vector to hold the observed measurements from all subjects
%indiv= []; % 1 x N vector to hold the indices of n subjects, e.g., tt(ind == 1), yy(ind == 1) for subject 1


%Initial out1 is based on the unique time points of the pooled data + the unique
%time points of "newdata", the output time grid. When newdata = [], output
%"out1" is equivalent to be the unique sorted pooled time points; otherwise, it 
%corresponds to the unique "newdata".

out1 =unique(tt);
out21 = linspace(min(out1),max(out1),ngrid);

mu = zeros(1,length(out1));

if strcmp(verbose, 'on') == 1
   fprintf(1,'Part II: Choose bandwidth of smoothing covariance surface\n');
end

rcov = getRawCov(y,t,out1,mu, regular, 0);           %obtain raw covariance;

if ~isempty(xcov)&&(length(xcov)==ngrid^2)
    if length(bwxcov)>1 && bwxcov(1)>0 && bwxcov(2)>0
        bw_xcov=bwxcov;
    else
        bw_xcov = [bw_mu bw_mu];
    end
else    
    if bwxcov(1)==0 || bwxcov(2)==0
        if (bwxcov_gcv == 1) || (bwxcov_gcv == 2)
            [bw_xcov] = gcv_mullwlsn(t,ngrid,regular,error, kernel,rcov,verbose);
            if isempty(bw_xcov)
                fprintf(1,'Error: FPCA is aborted because the observed data is too sparse to estimate the covariance function!');               
                return;
            end         
            bw_xcov = adjustBW2(kernel,bw_xcov,1,0,regular,verbose);
            if bwxcov_gcv == 2  
               minbwcov = getMinb(t,unique(cell2mat(t)),regular); %the minimum bw for cov function 
               bw_xcov = sqrt(minbwcov*bw_xcov);  
            end
            if strcmp(verbose, 'on') == 1
                 fprintf(1,['Adjusted GCV bandwidth choice for COV function : ' num2str(bw_xcov(1)) ',' num2str(bw_xcov(2)) ')\n']);
            end
        else
            [bw_xcov] = cv_mullwlsn(y,t,mu,ntest1,ngrid,regular,error,kernel,rcov,verbose);
        end
    elseif bwxcov > 0
        bw_xcov = bwxcov;
    elseif bwxcov(1) < 0 || bwxcov(2) < 0
        fprintf(1,'Error: Bandwidth choice for the covariance function must be positive!\n');
        return;
    end
    
    
    if isempty(xcov)
        rcov1 = rcov;
        if error == 1
            tpairn = rcov1.tpairn;
            tneq=find(tpairn(1,:)~=tpairn(2,:));
            cyy = rcov1.cyy;
            rcov1.tpairn = tpairn(:,tneq);
            rcov1.cxxn=cyy(tneq);
            rcov1.win=ones(1,length(rcov1.cxxn));
            if regular == 1
                rcov1.count = rcov1.count(tneq);
            end
        end
        
        if regular == 1
            [invalid,xcov]=mullwlsk(bw_xcov,kernel,rcov1.tpairn,rcov1.cxxn',rcov1.win,out21,out21,rcov1.count);  %smooth raw covariance;
        else
            [invalid,xcov]=mullwlsk(bw_xcov,kernel,rcov1.tpairn,rcov1.cxxn',rcov1.win,out21,out21);  %smooth raw covariance;
        end
    end
end 

xcov = (xcov+xcov')/2;   %transform the smoothed covariance matrix to guarantee it is a symmetric matrix.


end

function [regular] = isregular(t)
    tt = cell2mat(t);
    f = length(tt)/length(unique(tt))/length(t);
    if f==1
        regular = 2;
    elseif f>.75
        regular = 1;
    else
        regular = 0;
    end
end


