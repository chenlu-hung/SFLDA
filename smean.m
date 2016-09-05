
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
%Input bwmu:       scalar bwmu>=0, bandwidth for mean curve mu,
%                  0:  use cross-validation or generalized cross-validation
%                      to choose bandwidth automatically.      [Default]
%                  bwmu(>0): user-specified bandwidth.
%
%Input bwmu_gcv:   For choice bwmu = 0, two types of cross-validation can
%                  be performed: One-curve-leave-out cross-validation (CV)
%                  or  generalized cross-validation (GCV)
%                  0: CV   (may be time-consuming)
%                  1: Geometric mean of the minimum bandwidth and the GCV
%                     bandwidth.                            [Default]
%                  2: GCV  (faster than CV)
%

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




function [muDense out21 bw_mu]... 
         =smean(y,t,bwmu,bwmu_gcv,ntest1,ngrid,regular,kernel,numBins,verbose)

     no_opt =[];mu=[];muDense=[];bw_mu=[];out1 = [];out21 =[];xmu=[];
     

invalid = 0;
userrange = [];
[invalid_data]  = CheckData(y, t); %NaN values and number of subjects.
if (invalid_data == 1)
    return
end


%Set default values for the major input arguments
%when the following input arguments are set to 
%be "[]" by the user

if isempty(bwmu)
    bwmu = 0;    %bandwidth choice for mean function is using CV or GCV
end
if isempty(bwmu_gcv)
    bwmu_gcv = 1; % default bandwidth choice for mean function is the geometric mean of the minimum bandwidth and the GCV bandwidth
end

if isempty(ntest1)
    ntest1 = 30;
end
if isempty(ngrid)
   ngrid = 30;
end


if isempty(verbose)
  verbose = 'off';
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


if strcmp(verbose, 'on') == 1
  fprintf(1,'Part I: Obtain smoothed mean curve\n');
end

%when bwmu = 0 and bwmu_gcv = 0, use leave-one-curve-out CV method for bw choice
%when bwmu = 0 and bwmu_gcv = 1, use the geometric mean of the minimum bw and the leave-one-out GCV bw
%when bwmu = 0 and bwmu_gcv = 2, use leave-one-out GCV method for bw choice
%when bwmu > 0, user-defined bw for mean function

if ~isempty(xmu) && (length(xmu) == length(out1))
    mu = xmu;
    muDense = interp1(out1, mu, out21, 'spline');
    bw_mu = [];
else    
    if bwmu == 0
        if (bwmu_gcv == 1) || (bwmu_gcv == 2)
            bw_mu = gcv_lwls(yy,tt,kernel,1,1,0,regular,verbose);   %use GCV method to choose bw for mean function
            if isempty(bw_mu)
                fprintf(1,'Error: FPCA is aborted because the observed data is too sparse to estimate the mean function!');
                return;
            end
            bw_mu = adjustBW1(kernel,bw_mu,1,0,regular,verbose);
            if bwmu_gcv == 1  
                minbw = minb(tt,2); %the minimum bw for mean function 
                bw_mu = sqrt(minbw*bw_mu); %the geometric mean between the minimum bw and GCV bw
            end
            if strcmp(verbose, 'on') == 1
               fprintf(1,['Adjusted GCV bandwidth choice for mean function: ' num2str(bw_mu) '\n']);
            end
           
        else                                               %use CV method to choose bw for mean function
            bw_mu = cvfda_lwls(y,t,kernel,1,1,0,regular,verbose);
        end
    elseif bwmu > 0
        bw_mu = bwmu;
    else
        fprintf(1,'Error: Bandwidth choice for the mean function must be positive!\n');
        return;
    end
    
    %define the vector of case weight in the local weighted least square
    %here, it is set to be one for all subjects
    % win1 = ones(1,length(tt));
    win = cell(1,length(t));
    for i = 1:length(t)
        win{i} = ones(1,ni(i))/ni(i);
    end
    win1 = cell2mat(win);
    [invalid, mu] = lwls(bw_mu,kernel,1,1,0,tt,yy',win1,out1);
    [invalid, muDense] = lwls(bw_mu,kernel,1,1,0,tt,yy',win1,out21);
    
end



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


