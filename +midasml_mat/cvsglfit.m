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
% Updated: 20240213

p = inputParser;
addRequired(p,'x', @(z) isnumeric(z));
addRequired(p,'y', @(z) isnumeric(z));
addParameter(p,'gamma',1, @(z) isnumeric(z));
addParameter(p,'nlambda',100, @(z) isnumeric(z));   
addParameter(p,'lambda_factor',[], @(z) isnumeric(z)); 
addParameter(p,'lambda',[], @(z) isnumeric(z));           
addParameter(p,'pf',[], @(z) isnumeric(z)); 
addParameter(p,'gindex',[], @(z) isnumeric(z)) 
addParameter(p,'dfmax',[], @isnumeric);        
addParameter(p,'pmax',[], @(z) isnumeric(z)); 
addParameter(p,'standardize',false, @(z) islogical(z)); 
addParameter(p,'intercept',true, @(z) islogical(z)); 
addParameter(p,'eps',1e-8, @(z) isnumeric(z));  
addParameter(p,'maxit',1e6, @(z) isnumeric(z)); 
addParameter(p,'peps',1e-8, @(z) isnumeric(z)); 
addParameter(p,'fe',false, @(z) islogical(z)); 
addParameter(p,'N',[], @(z) isnumeric(z)); 
addParameter(p,'nfolds',10, @(z) isnumeric(z)); 
addParameter(p,'foldid',[], @(z) isnumeric(z)); 
addParameter(p,'parallel',false, @(z) islogical(z)); 

parse(p,x,y,varargin{:});
y = p.Results.y;
x = p.Results.x;
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
nfolds = p.Results.nfolds; 
foldid = p.Results.foldid;
parallel = p.Results.parallel;

% [nfolds,foldid,parallel,gamma,nlambda,lambda_factor,lambda,pf,gindex,dfmax,pmax,standardize,intercept,eps,maxit,peps,fe,N] = ...
%     process_options(varargin,'nfolds',10,'foldid',[],'parallel',false,'gamma',1.0,'nlambda',100,'lambda_factor',[],'lambda',[],'pf',[],'gindex',[], ...
%     'dfmax',[],'pmax',[],'standardize',true,'intercept',true,'eps',1e-8,'maxit',1e6,'peps',1e-8,'fe',false,'N',[]);

if isempty(nfolds)
    nfolds=10;
end

[n,~] = size(x);

if isempty(foldid)
    population = cat(2, repmat(1:nfolds, 1, floor(n/nfolds)), 1:mod(n,nfolds));
    foldid = sort(population);
end
nfolds = max(foldid);

sglfit = sgl(x,y,'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
    'lambda',lambda,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
    'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',fe,'N',[]);
lambdaz = sglfit.lambda; %The lambda grid

cvraw = zeros(nfolds, nlambda);
if fe %PANEL (fe)
    if parallel==true %PARALLEL
        parfor i = 1:nfolds %PARFOR-loop
            which = (foldid == i);
            xin = x(~which,:);
            yin = y(~which,:);
            fit =  sgl(xin,yin,'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
                'lambda',lambdaz,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
                'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',true,'N',N); %here: fe=true
            e_lam = sglfitpredict(fit,  x(which,:),fe) - repmat(y(which,:), 1, nlambda); 
            cvraw(i, :) = mean(e_lam.^2,1);
        end
    else %SERIAL (NOT-PARALLEL)
        for i = 1:nfolds %FOR-loop
            which = (foldid == i);
            fit =  sgl(x(~which,:),y(~which,:),'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
                'lambda',lambdaz,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
                'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',true,'N',N); %here: fe=true
            e_lam = sglfitpredict(fit,  x(which,:),fe) - repmat(y(which,:), 1, nlambda);
            cvraw(i, :) = mean(e_lam.^2,1);
        end
    end
    
else %NOT PANEL (fe=false)
    if parallel==true
        parfor i = 1:nfolds
            which = (foldid == i);
            xin = x(~which,:);
            yin = y(~which,:);
            fit =  sgl(xin,yin,'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
                'lambda',lambdaz,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
                'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',false,'N',[]); %here: fe=false
            e_lam = sglfitpredict(fit,  x(which,:),fe) - repmat(y(which,:), 1, nlambda);
            cvraw(i, :) = mean(e_lam.^2,1);
        end
    else
        for i = 1:nfolds
            %disp(['Currently doing fold i=: ' num2str(i) ' out of ' num2str(nfolds)])
            which = (foldid == i);
            fit =  sgl(x(~which,:),y(~which,:),'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
                'lambda',lambdaz,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
                'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',false,'N',[]); %here: fe=false
            e_lam = sglfitpredict(fit,  x(which,:),fe) - repmat(y(which,:), 1, nlambda);
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
%lambda_min
output.lambda_min = max(sglfit.lambda(cvm<=min(cvm)));
%lambda_1se
idmin = (sglfit.lambda==output.lambda_min);
semin = cvm(idmin)+cvsd(idmin);
output.lambda_1se = max(sglfit.lambda(cvm<=semin));
idmin_1se = (sglfit.lambda==output.lambda_1se);
%betas
output.cvsglfit.lam_min.b0 = sglfit.b0(idmin);
output.cvsglfit.lam_min.beta = sglfit.beta(:,idmin);
output.cvsglfit.lam_1se.b0 = sglfit.b0(idmin_1se);
output.cvsglfit.lam_1se.beta = sglfit.beta(:,idmin_1se);
output.class = 'cv.sglfit';
end