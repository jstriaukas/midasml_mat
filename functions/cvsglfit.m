function output = cvsglfit(x,y,varargin)
%--------------------------------------------------------------------------
% cvsglfit: fit a linear model with sg-LASSO regularization
%           and return cross-validation solution.
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
% 'nfolds'         number of folds of the cv loop. Default set to 10.
% 'foldid'         the fold assignments used.
% 'parallel'       should parfor be used in a cv loop. Default set to false
% 'gamma'          sg-LASSO mixing parameter. γ = 1 gives LASSO solution and γ = 0 gives group LASSO solution.
%
%  FOr other optional arguments, see sgl.m
%
% OUTPUT ARGUMENTS:
%   output     List of elements of cvsglfit.
%   output.sgl_fit     Full sample sg-LASSO fit.
%   output.cvm         Average of CV error curve.
%   output.cvsd        Standard deviation of cross validation error curve.
%   output.lambda_min  lambda which minimizes CV error curve.
%   output.lambda_1se  lambda_min + 1se.
%   output.cvsglfit    CV fit.
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
[nfolds,foldid,parallel,gamma,nlambda,lambda_factor,lambda,pf,gindex,dfmax,pmax,standardize,intercept,eps,maxit,peps,fe,N] = ...
    process_options(varargin,'nfolds',10,'foldid',[],'parallel',false,'gamma',1.0,'nlambda',100,'lambda_factor',[],'lambda',[],'pf',[],'gindex',[], ...
    'dfmax',[],'pmax',[],'standardize',false,'intercept',true,'eps',1e-8,'maxit',1e6,'peps',1e-8,'fe',false,'N',[]);

[n,~] = size(x);

if isempty(foldid)
    population = cat(2, repmat(1:nfolds, 1, floor(n/nfolds)), 1:mod(n,nfolds));
    foldid = sort(population);
end

sglfit = sgl(x,y,'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
    'lambda',lambda,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
    'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps);
cvraw = zeros(nfolds, nlambda);
if fe
    if parallel==true
        parfor i = 1:nfolds
            which = foldid == i;
            xin = x(~which,:);
            yin = y(~which,:);
            fit =  sgl(xin,yin,'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
                'lambda',lambda,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
                'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',true,'N',N);
            e_lam = sglfitpredict(fit,  x(which,:)) - repmat(y(which,:), 1, nlambda);
            cvraw(i, :) = mean(e_lam.^2,1);
        end
    else
        for i = 1:nfolds
            which = foldid == i;
            fit =  sgl(x(~which,:),y(~which,:),'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
                'lambda',lambda,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
                'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',true,'N',N);
            e_lam = sglfitpredict(fit,  x(which,:)) - repmat(y(which,:), 1, nlambda);
            cvraw(i, :) = mean(e_lam.^2,1);
        end
    end
    
else
    
    if parallel==true
        parfor i = 1:nfolds
            which = foldid == i;
            xin = x(~which,:);
            yin = y(~which,:);
            fit =  sgl(xin,yin,'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
                'lambda',lambda,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
                'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps);
            e_lam = sglfitpredict(fit,  x(which,:)) - repmat(y(which,:), 1, nlambda);
            cvraw(i, :) = mean(e_lam.^2,1);
        end
    else
        for i = 1:nfolds
            which = foldid == i;
            fit =  sgl(x(~which,:),y(~which,:),'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
                'lambda',lambda,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
                'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps);
            e_lam = sglfitpredict(fit,  x(which,:)) - repmat(y(which,:), 1, nlambda);
            cvraw(i, :) = mean(e_lam.^2,1);
        end
    end
    
end

cvm = mean(cvraw,1);
sqccv = (bsxfun(@minus,cvraw,cvm)).^2;
cvsd = sqrt(mean(sqccv,1)./(nfolds-1));
output.sglfit = sglfit;
output.cvm = cvm';
output.cvsd = cvsd';
output.lambda_min = max(sglfit.lambda(cvm<=min(cvm)));
idmin = sglfit.lambda==output.lambda_min;
semin = cvm(idmin)+cvsd(idmin);
output.lambda_1se = max(sglfit.lambda(cvm<=semin));
output.cvsglfit.lam_min.b0 = sglfit.b0(idmin);
output.cvsglfit.lam_min.beta = sglfit.beta(:,idmin);
idmin_1se = sglfit.lambda==output.lambda_1se;
output.cvsglfit.lam_1se.b0 = sglfit.b0(idmin_1se);
output.cvsglfit.lam_1se.beta = sglfit.beta(:,idmin_1se);
output.class = 'cv.sglfit';
end