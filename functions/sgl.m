function output = sgl(x,y,varargin)
%--------------------------------------------------------------------------
% sgl: fit a linear model with sg-LASSO regularization.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Fit a (currently only) linear model via penalized least squares. The 
%    regularization path is computed for the sg-LASSO at a grid of values 
%    for the regularization parameters lambda and gamma. 
%
% USAGE:
%   output = sgl(x,y)
%
% INPUT ARGUMENTS:
% x           Input matrix, of dimension nobs x nvars; each row is an
%             observation vector.
% y           Response variable. 
%
% OPTIONAL ARGUMENTS: 
% 'gamma'          sg-LASSO mixing parameter. γ = 1 gives LASSO solution and γ = 0 gives group LASSO solution.
% 'nlambda'        number of λ’s to use in the regularization path; used if lambda = NULL
% 'lambda_factor'  The factor for getting the minimal λ in the λ sequence, where min(lambda) =
%                  lambda.factor * max(lambda). max(lambda) is the smallest value of lambda
%                  for which all coefficients are zero. λmax is determined for each γ tuning parameter separately. 
%                  The default depends on the relationship between T (the sample size) and p (the number of 
%                  predictors). If T < p, the default is 0.01. If T > p, the default is 0.0001, closer to zero. 
%                  The smaller the value of lambda.factor is, the denser is the fit for λmin. Used only if lambda = NULL.
% 'lambda'         a user-supplied lambda sequence. By leaving this option unspecified (recommended), users can have the 
%                  program compute its own lambda sequence based on nlambda and lambda.factor. It is better to supply, 
%                  if necessary, a decreasing sequence of lambda values than a single (small) value, as warm-starts are
%                  used in the optimization algorithm. The program will ensure that the usersupplied λ sequence is sorted 
%                  in decreasing order before fitting the model.
% 'pf'             the ell_1 penalty factor of length p used for the adaptive sg-LASSO. Separate ell_1 penalty weights 
%                  can be applied to each coefficient to allow different ell_1 + ell_2,1 shrinkage. Can be 0 for some 
%                  variables, which imposes no shrinkage, and results in that variable always be included in the model. 
%                  Default is 1 for all variables.
% 'gindex'         p by 1 vector indicating group membership of each covariate.              
% 'dfmax'          the maximum number of variables allowed in the model. Useful for very large
%                  p when a partial path is desired. Default is p+1. In case method='fe', dfmax is ignored.
% 'pmax'           the maximum number of coefficients allowed ever to be nonzero. For example, once βi ~= 0 for some i ∈ [p], 
%                  no matter how many times it exits or re-enters the model through the path, it will be counted only once. 
%                  Default is min(dfmax*1.2,p).
% 'standardize'    logical flag for variable standardization, prior to fitting the model sequence. The
%                  coefficients are always returned to the original scale. It is recommended to keep
%                  standardize=TRUE. Default is FALSE.     
% 'intercept'      whether intercept be fitted (TRUE) or set to zero (FALSE). Default is FALSE.                
% 'eps'            convergence threshold for block coordinate descent. Each inner block coordinatedescent loop continues until
%                  the maximum change in the objective after any coefficient update is less than thresh times the null deviance. 
%                  Defaults value is 1e-8.
% 'maxit'          maximum number of outer-loop iterations allowed at fixed lambda values. Default is 1e6. 
%                  If the algorithm does not converge, consider increasing maxit.
% 'peps'           convergence threshold for proximal map of sg-LASSO penalty. Each loop continues until G group difference 
%                  sup-norm, ||β^k_G - β^k-1_G||_∞, is less than peps. Defaults value is 1e-8     
% 'fe'             If we need to fit fixed effects. 'N', number of fixed
%                  effects must be supplied in this case.
% 'N'              Number of fixed effects in panel regression.
%          
%
% DETAILS:
%              Estimates sg-LASSO model under MSE loss:
%            
%              min_b |y-Xb|_T^2 + 2lambda (gamma |b|_1 + (1-gamma)|b|_2,1.  
%
%              A pair of tuning parameters (lambda,gamma) is chosen via
%              10-fold cross-validation. 
%
% DATE: 2021-06-09
% AUTHOR: Jonas Striaukas
% LICENSE: GPL-2
[gamma,nlambda,lambda_factor,lambda,pf,gindex,dfmax,pmax,standardize,intercept,eps,maxit,peps,fe,N] = ...
        process_options(varargin,'gamma',1,'nlambda',100,'lambda_factor',[],'lambda',[],'pf',[],'gindex',[], ...
        'dfmax',[],'pmax',[],'standardize',false,'intercept',true,'eps',1e-8,'maxit',1e6,'peps',1e-8,'fe',false,'N',[]);

nobs = int32(size(x,1));
nvars = int32(size(x,2));
if (isempty(lambda_factor))
    if (nobs < nvars)
        lambda_factor = 1e-2;
    else
        lambda_factor = 1e-4;
    end
end
maxit = int32(maxit);
nlam = int32(nlambda);
if (isempty(gindex))
    gindex = 1:nvars;
end
if (size(gindex,1)==1)
    gindex = gindex';
end
if (isempty(pf))
    pf = ones(nvars, 1);
end
if (size(pf,1)==1)
    pf = pf';
end
if (isempty(dfmax))
    dfmax = nvars+1;
end
if (isempty(pmax))
    pmax = min(1.2 * dfmax, nvars);
end

ngroups = int32(max(gindex));
if (any(diff(gindex)>1))
    error('only adjancet group memberships allowed');
end

if (isempty(lambda))
    if (lambda_factor >= 1.0)
        error('lambda_factor should be less than 1.0')
    end
    flmin = double(lambda_factor);
    ulam = computelambda(nlambda, flmin, nobs, x, y, gamma, gindex, length(gindex), pf);
else
    flmin = double(1);
    ulam = double(sort(lambda, 'descend'));
    nlam = int32(length(lambda));
end
gamma = double(gamma);
if (gamma < 0 || gamma > 1)
    error('gamma must be in [0,1]');
end

gindex = int32(find([diff(gindex);1]==1));
nobs = int32(nobs);
nvars = int32(nvars);
x = double(x);
y = double(y);
pf = double(pf);
dfmax = int32(dfmax);
pmax = int32(pmax);
nlam = int32(nlam);
eps = double(eps);
peps = double(peps);
isd = int32(standardize);
intr = int32(intercept);
maxit = int32(maxit);

load linux_output.mat
nalam1 = nalam;
b01 = b0;
beta1 = beta;
ibeta1 = ibeta;
nbeta1 = nbeta;
alam1 = alam;
npass1 = npass;
jerr1 = jerr;
% --------------------------------- fit sg-LASSO -------------------------%
if fe
    yn = y;
    xn = x;
    T = size(x,1)/N;
    if rem(T,1)~=0
        error('you chose to fit fixed effects with sg-LASSO, but the number of fixed effects specified is wrong, i.e. NT != N * T. change ''N'' or set ''fe=false''.')
    end
    ymb = zeros(N,1);
    xmb = zeros(N,size(x,2));
    for k = 1:N
        festart = (k-1)*T+1;
        feend = k*T;
        ymb(k) = mean(y(festart:feend));
        xmb(k,:) = mean(x(festart:feend,:),1);
        yn(festart:feend) = y(festart:feend) - ymb(k);
        xn(festart:feend,:) = x(festart:feend,:) - xmb(k,:);
    end
    intr = int32(0);
     [nalam,b0,beta,ibeta,nbeta,alam,npass,jerr] = sglfit(gamma, ngroups, gindex, ...
            nobs, nvars, xn, yn, pf, dfmax, pmax, nlam, flmin, ulam, ... 
            eps, peps, isd, int32(intr), maxit);
else
    
    [nalam,b0,beta,ibeta,nbeta,alam,npass,jerr] = sglfit(gamma, ngroups, gindex, ...
            nobs, nvars, x, y, pf, dfmax, pmax, nlam, flmin, ulam, ... 
            eps, peps, isd, int32(intr), maxit);
end
% ------------------------------------------------------------------------%   
if (jerr~=0)
    errmsg = err(jerr,maxit,pmax);
    if (errmsg.fatal)
        error(errmsg.msg);
    else
        warning(errmsg.msg);
    end
end
b = zeros(nvars, nlam);
for l = 1:nalam
    nk = nbeta(l);
    b(ibeta(1:nk),l) = beta(1:nk,l);
end
if fe
    b0 = zeros(nlambda, 1);
    a0 = zeros(N, nlambda);
    for l = 1:nlambda
        a0(:,l) = ymb - xmb*b(:,l);
    end
    output.a0 = a0;
    output.N = N;
end
output.b0 = b0;
output.beta = b;
output.npass = npass;
output.lambda = alam;

end

% ----------------- error handling functions -----------------------%
function output = err(n,maxit,pmax)

    if n==0
        output.n=0;
        output.fatal=false;
        output.msg='';
    else
        output = errsgl(n,maxit,pmax);    
        output.msg = sprintf('from glmnet Fortran code (error code %d); %s', n, output.msg);
    end
end

function output = errsgl(n,maxit,pmax)

    if (n > 0)  %fatal error
        if (n < 7777)
            msg = 'Memory allocation error; contact package maintainer';
        elseif (n == 7777)
            msg = 'All used predictors have zero variance';
        elseif (n == 10000)
            msg = 'All penalty factors are <= 0';
        else
            msg = 'Unknown error';
        end
        output.n = n;
        output.fatal = true;
        output.msg = msg;
    elseif (n < 0)  %non-fatal error
        if (n > -10000)
            msg = sprintf('Convergence for %dth lambda value not reached after maxit=%d iterations; solutions for larger lambdas returned',-n,maxit);
        end
        if (n < -10000)
            msg = sprintf('Number of nonzero coefficients along the path exceeds pmax=%d at %dth lambda value; solutions for larger lambdas returned',pmax,-n-10000);
        end
        output.n = n;
        output.fatal = false;
        output.msg = msg;
    end
end


    
    