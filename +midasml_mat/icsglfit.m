function output = icsglfit(x,y,ic,varargin)
%--------------------------------------------------------------------------
% cvsglfit: fit a linear model with sg-LASSO regularization 
%           and return solution based on information criteria (IC) 
%           choice.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Fit a (currently only) linear model via penalized least squares. The 
%    regularization path is computed for the sg-LASSO at a grid of values 
%    for the regularization parameters lambda and gamma. 
%
% USAGE:
%   output = cvsglfit(x,y)
%
% INPUT ARGUMENTS:
% x           Input matrix, of dimension nobs x nvars; each row is an
%             observation vector.
% y           Response variable. 
%
% OPTIONAL ARGUMENTS: 
% 'gindex'         p by 1 vector indicating group membership of each covariate.
% 'ic'             choice of information criteria (AIC (aic), BIC (bic) 
%                   or AICc (aicc)). Default BIC.
% 'gamma'          sg-LASSO mixing parameter. γ = 1 gives LASSO solution and γ = 0 gives group LASSO solution.
% 
%  FOr other optional arguments, see sgl.m
%
% OUTPUT ARGUMENTS:
%   output     List of elements of cvsglfit.
%   output.sgl_fit     Full sample sg-LASSO fit.
%   output.cvm         Average of CV error curve.
%   output.lambda_ic  lambda which minimizes IC criteria.
%   output.icsglfit   IC fit.
%
% DETAILS:
%              Estimates sg-LASSO model under MSE loss:
%            
%              min_b |y-Xb|_T^2 + 2lambda (gamma |b|_1 + (1-gamma)|b|_2,1.  
%
%              A pair of tuning parameters (lambda,gamma) is chosen via
%              10-fold cross-validation. 
%
% DATE: 2021-11-10
% AUTHOR: Jonas Striaukas
% LICENSE: GPL-2
% Updated: 20240213

p = inputParser;
addRequired(p,'x', @(z) isnumeric(z));
addRequired(p,'y', @(z) isnumeric(z));
addRequired(p,'ic',@(z) islogical(z) || ischar(z) || isnumeric(z)); 
addParameter(p,'gamma',1, @(z) isnumeric(z)); 
addParameter(p,'nlambda',100, @(z) isnumeric(z));         
addParameter(p,'lambda_factor',[], @(z) isnumeric(z)); 
addParameter(p,'lambda',[], @(z) isnumeric(z));           
addParameter(p,'pf',[], @(z) isnumeric(z)); 
addParameter(p,'gindex',[], @(z) isnumeric(z)) 
addParameter(p,'dfmax',[], @(z) isnumeric(z));     
addParameter(p,'pmax',[], @(z) isnumeric(z)); 
addParameter(p,'standardize',false, @(z) islogical(z)); 
addParameter(p,'intercept',true, @(z) islogical(z)); 
addParameter(p,'eps',1e-8, @(z) isnumeric(z));
addParameter(p,'maxit',1e6, @(z) isnumeric(z));
addParameter(p,'peps',1e-8, @(z) isnumeric(z));
addParameter(p,'fe',false, @(z) islogical(z)); 
addParameter(p,'N',[], @(z) isnumeric(z)); 

parse(p,x,y,ic,varargin{:});
x = p.Results.x;
y = p.Results.y;
ic = p.Results.ic;
gamma = p.Results.gamma;
nlambda = p.Results.nlambda;
lambda_factor = p.Results.lambda_factor;
lambda = p.Results.lambda;
pf = p.Results.pf;
gindex = p.Results.gindex;
dfmax = p.Results.dfmax;
pmax = p.Results.pmax;
standardize = p.Results.standardize; 
intercept = p.Results.intercept;
eps = p.Results.eps;
maxit = p.Results.maxit;
peps = p.Results.peps;
fe = p.Results.fe;
N = p.Results.N;

%     [ic,gamma,nlambda,lambda_factor,lambda,pf,gindex,dfmax,pmax,standardize,intercept,eps,maxit,peps] = ...
%             process_options(varargin,'ic','bic','gamma',1.0,'nlambda',100,'lambda_factor',[],'lambda',[],'pf',[],'gindex',[], ...
%             'dfmax',[],'pmax',[],'standardize',false,'intercept',true,'eps',1e-8,'maxit',1e6,'peps',1e-8);

    ic = char(lower(ic)); %string containing the ic: e.g. 'BIC'

    if isempty(ic) %Set default ic value
        ic='bic'; 
    end
    
    validICs = {'bic','aic','aicc'};
    if ~any(strcmpi(ic, validICs))
        error(['Input ''ic'' can only be one of the following: ', char(join(validICs,", "))] )
    end
    
    [n,~] = size(x);
    sglfit = sgl(x,y,'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
                'lambda',lambda,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
                'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',false,'N',N); %here: fe=false
    
    fits = sglfitpredict(sglfit, x); %the fitted_y's for all the lambdas
    cvm = zeros(nlambda,1); %preallocation
    sigsqhat = sum((y-mean(y)).^2)/n;
    for i = 1:nlambda %Loop over lambdas
        mse = mean((y - fits(:,i)).^2)/n;
        df = sum(sglfit.beta(:,i) == 0) + double(intercept);
        cvm(i) = mse/sigsqhat + ic_pen(ic, df, n); %BIC, AIC, or AICc for each lambda
    end
    
    output.sglfit = sglfit;
    output.cvm = cvm'; %All the ICs
    idxmin = find(cvm==min(cvm),1); %NB: If the smallest element occurs multiple times, takes the 1st one (and since lambda sequence is decsending, this takes the most parsimonious solution)

    output.lambda_min = sglfit.lambda(idxmin); %The optimal LAMBDA
    output.icsglfit.lam_min.b0 = sglfit.b0(idxmin);
    output.icsglfit.lam_min.beta = sglfit.beta(:,idxmin);

    output.class = 'cv.sglfit';
end

%%% Subfunction %%%
function pen = ic_pen(ic_choice, df, t)
  if strcmp(ic_choice,"bic") %BIC
      pen = log(t)/t*df;
  elseif strcmp(ic_choice,"aic") %AIC
      pen = 2/t*df;
  elseif strcmp(ic_choice,"aicc") %AICc
      pen = 2*df/(t - df - 1);
  end
end
