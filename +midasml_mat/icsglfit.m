function output = icsglfit(x,y,varargin)
%--------------------------------------------------------------------------
% icsglfit: fit a linear model with sg-LASSO regularization 
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
    [ic,gamma,nlambda,lambda_factor,lambda,pf,gindex,dfmax,pmax,standardize,intercept,eps,maxit,peps] = ...
            process_options(varargin,'ic','bic','gamma',1.0,'nlambda',100,'lambda_factor',[],'lambda',[],'pf',[],'gindex',[], ...
            'dfmax',[],'pmax',[],'standardize',false,'intercept',true,'eps',1e-8,'maxit',1e6,'peps',1e-8);
    [n,~] = size(x);
    sglfit = sgl(x,y,'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
                'lambda',lambda,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
                'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps);
    fits = sglfitpredict(sglfit, x);
    cvm = zeros(nlambda,1);
    sigsqhat = sum((y-mean(y)).^2)/n;
    for i = 1:nlambda
        mse = mean((y - fits(:,i)).^2)/n;
        df = sum(sglfit.beta(:,i) == 0) + double(intercept);
        cvm(i) = mse/sigsqhat + ic_pen(ic, df, n);
    end
    output.sglfit = sglfit;
    output.cvm = cvm';
    idxmin = min(cvm)==cvm;
    output.lambda_min = sglfit.lambda(idxmin);
    output.icsglfit.lam_min.b0 = sglfit.b0(idxmin);
    output.icsglfit.lam_min.beta = sglfit.beta(:,idxmin);

    output.class = 'ic.sglfit';
end

function pen = ic_pen(ic_choice, df, t)
  if strcmp(ic_choice,"bic")
    pen = log(t)/t*df;
  end
  if strcmp(ic_choice,"aic")
    pen = 2/t*df;
  end
  if strcmp(ic_choice,"aicc")
    pen = 2*df/(t - df - 1);
  end
end
