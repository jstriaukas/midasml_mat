function output = tscvsglfit(x,y,varargin)

%INPUTS:
%- 'l' is the gap used to drop observations around the test set data point. 
%  Default: l=5.
%- 'K' number of observations drawn for the test set. 
%  Default: K=20
%- X is T-by-Nx
%- For all other inputs see the sgl() function inputs.

%Example:
%obj = tscvsglfit(x,y,'gamma',0.5,'gindex',gindex,'K',20,'l',5);
%[yhat,bz] = cvsglfitpredict(obj, newX, 'lambda_1se')

% Updated: 20240213

p = inputParser;
addRequired(p,'x',@(z) isnumeric(z) || isstruct(z));
addRequired(p,'y',@isnumeric);
addParameter(p,'parallel',false,@(z) islogical(z)); 
addParameter(p,'gamma',1,@(z) isnumeric(z)); 
addParameter(p,'nlambda',100,@(z) isnumeric(z)); 
addParameter(p,'lambda_factor',[],@(z) isnumeric(z)); 
addParameter(p,'lambda',[],@(z) isnumeric(z)); 
addParameter(p,'pf',[],@(z) isnumeric(z)); 
addParameter(p,'gindex',[],@(z) isnumeric(z)); 
addParameter(p,'dfmax',[],@(z) isnumeric(z)); 
addParameter(p,'pmax',[],@(z) isnumeric(z)); 
addParameter(p,'standardize',true,@(z) islogical(z)); 
addParameter(p,'intercept',true,@(z) islogical(z)); 
addParameter(p,'eps',1e-8,@(z) isnumeric(z)); 
addParameter(p,'maxit',1e6,@(z) isnumeric(z)); 
addParameter(p,'peps',1e-8,@(z) isnumeric(z)); 
addParameter(p,'fe',false,@(z) islogical(z)); 
addParameter(p,'l',5,@(z) isnumeric(z)); 
addParameter(p,'K',20,@(z) isnumeric(z)); 

parse(p,x,y,varargin{:});
x = p.Results.x;
y = p.Results.y;
parallel = p.Results.parallel;
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
fe = p.Results.fe; %If set to ture, returns error
l = p.Results.l;
K = p.Results.K;

if fe
   error('Name-value pair ''fe'' cannot be set to true when using TS-CV, as it''s strictly for time series data.') 
end

if parallel
    parallel = false;
    warning('tscvsglfit is currnetly only running in serial. ''Parallel'' was set to false, automatically.')
end

[T,~] = size(x); 
if l<=1
   error('Input ''l'' must satisfy: l>=2') 
end
if K<1 || K>=T
    error('Input ''k'' must satisfy: 1<= K <= T-1')
end

sglfit = sgl(x,y,'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
    'lambda',lambda,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
    'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',fe);
lambdaz = sglfit.lambda; %The lambda grid

foldid = randsample(1:T, K, false); %sample K values from 1:T without replacement
o = cell(K,1); 
L = length(lambdaz);
yhat = nan(length(y), L);
cvraw = nan(K,L);

if ~parallel
    for i = 1:K 
        idxtest = foldid(i); %TEST set
        idxtrain = setdiff(1:T, (idxtest-l:idxtest+l)); %TRAINING set
        o{i} = sgl(x(idxtrain,:),y(idxtrain),'gamma',gamma,'nlambda',nlambda,'lambda_factor',lambda_factor,...
            'lambda',lambdaz,'pf',pf,'gindex',gindex,'dfmax',dfmax,'pmax',pmax,'standardize',standardize,...
            'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',fe);

        yhat(idxtest,:) = sglfitpredict(o{i},  x(idxtest,:), false); %T-by-L
        e_lam = yhat(idxtest,:) - repmat(y(idxtest), [1 L]); %(Ttest x L)
        cvraw(i,:) = mean(e_lam.^2,1); %(K x L)
    end
else 

end
cvm = mean(cvraw,1); %(1xL)
cvsd  = std(cvraw)/sqrt(size(cvraw,1));

output.sglfit = sglfit;
output.cvm = cvm';
output.cvsd = cvsd';
%lambda_min
output.lambda_min = max(sglfit.lambda(cvm<=min(cvm)));
idmin = (sglfit.lambda==output.lambda_min);
%lambda_1se
semin = cvm(idmin)+cvsd(idmin);
output.lambda_1se = max(sglfit.lambda(cvm<=semin));
idmin_1se = (sglfit.lambda==output.lambda_1se);
%betas
output.cvsglfit.lam_min.b0 = sglfit.b0(idmin);
output.cvsglfit.lam_min.beta = sglfit.beta(:,idmin);
output.cvsglfit.lam_1se.b0 = sglfit.b0(idmin_1se);
output.cvsglfit.lam_1se.beta = sglfit.beta(:,idmin_1se);
output.class = 'tscv.sglfit';

end