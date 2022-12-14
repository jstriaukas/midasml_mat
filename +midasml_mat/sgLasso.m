function probj = sgLasso(Y,X,varargin)
    
    %sgLasso() estimates Sparse Group LASSO (sg-LASSO) regression. It allows 
    %to select the optimal tunning parameters Lambda & Gamma using either 
    %information criteria (BIC, AIC, and AICc) or Kfold-CV.
    %- How Hyperparameter optimization works in sgLasso():
    %  - When Gamma is a single value (gamma=0.5) -> HPOptimizeOptions.Method is used to pick the optimal Lambda
    %  - When Gamma is a grid (gamma=0:0.1:1)     -> HPOptimizeOptions.Method is used to pick the optimal Lambda & Gamma
    %
    %INPUTS:
    %- 'Y' is a Tx1 vector
    %- 'X' is a TxN matrix. X should NOT contain vector of ones for the
    %  constant. A constant is added automatically internally.
    %- 'Gamma' is the LASSO/Group-LASSO mixing parameter, and it can be either: 
    %  1. a SCALAR: This will estimate the model at that fixed value of gamma.
    %  2. a VECTOR: This will optimize over the grid of provided gamma values.
    %  The two extremes correspond to Group-LASSO (gamma=0) and LASSO (gamma=1).
    %
    %OUTPUT:
    %fit.b <- (1+K)x1 vector with the betas (including constant in 1st row)
    %          where K=N+Nd+Ny; N=size(X,2); Nd=size(D,2); Ny=size(Ylags)
    
    %Examples:
    %
    %Some data
    %T=1000; N=20;
    %X = [ones(T,1) randn(T,N)]; 
    %b = zeros(N+1,1); 
    %b(2:6) = 2;     
    %b(7:11) = -3;
    %D=[];
    %y = X*b + randn(T,1)*0.1;
    %x = X(:,2:end); 
    %clearvars -except y x
    %
    %o1 = struct('Method',{'BIC'})
    %o2 = struct('Method',{'AIC'})
    %o3 = struct('Method',{'CV-Kfold'}, 'Folds',5, 'OptimMeth',{'MinErr'})
    %o4 = struct('Method',{'CV-Kfold'}, 'Folds',10, 'OptimMeth',{'MinErrPlus1SE'})
    %gindex = [];
    %
    %- gamma=0.5 and optimal Lambda selected via CV-Kfold
    %fit = sgLasso(y,x,'ctrl',[],'Groupctrl',false,'ARlags',[],'GroupARlags',false,'Gamma',0.5,'gindex',gindex,'HPOptimizeOptions',o3)
    %
    %- Optimize both Gamma & Lambda via BIC
    %gamma = linspace(1,0,21);
    %fit = sgLasso(y,x,'ctrl',[],'Groupctrl',false,'ARlags',[],'GroupARlags',false,'Gamma',gamma,'gindex',gindex,'HPOptimizeOptions',o1)
    %fit.b %<- The estimates betas (inc. the constant)

    p = inputParser;
    addRequired(p,'Y',@isnumeric);
    addRequired(p,'X',@(z) isnumeric(z) || isstruct(z));
    addParameter(p,'Gamma',1,@(z) isnumeric(z)); 
    addParameter(p,'ctrl',[],@isnumeric)
    defoo = struct('Method',{'CV-Kfold'}, 'Folds',10, 'OptimMeth',{'MinErr'});
    addParameter(p,'HPOptimizeOptions',defoo,@(z) isstruct(z));
    addParameter(p,'parallel',false,@islogical)
    addParameter(p,'standardize',false,@(z) islogical(z));
    addParameter(p,'dfmax',[],@isnumeric)
    addParameter(p,'ARlags',[],@isnumeric); 
    addParameter(p,'GroupARlags',false,@(z) islogical(z));
    addParameter(p,'Groupctrl',false,@(z) islogical(z));
    addParameter(p,'gindex',[],@(z) isnumeric(z))
    %addParameter(p,'nfolds',10,@(z) isnumeric(z)); 
    %addParameter(p,'foldid',[],@(z) isnumeric(z)); 
    %addParameter(p,'nlambda',100,@(z) isnumeric(z)); 
    %addParameter(p,'lambda_factor',[],@(z) isnumeric(z)); 
    %addParameter(p,'lambda',[],@(z) isnumeric(z)); 
    addParameter(p,'pf',[],@(z) isnumeric(z)); 
    addParameter(p,'pmax',[],@(z) isnumeric(z)); 
    addParameter(p,'intercept',true,@(z) islogical(z)); 
    addParameter(p,'eps',1e-8,@(z) isnumeric(z)); 
    addParameter(p,'maxit',1e6,@(z) isnumeric(z)); 
    addParameter(p,'peps',1e-8,@(z) isnumeric(z)); 
    addParameter(p,'fe',false,@(z) islogical(z)); 
    addParameter(p,'N',[],@(z) isnumeric(z)); 

    parse(p,Y,X,varargin{:});
    Y = p.Results.Y;
    X = p.Results.X;
    D = p.Results.ctrl;
    alfa = p.Results.Gamma;
    oo = p.Results.HPOptimizeOptions;
    parallel = p.Results.parallel;
    standardize = p.Results.standardize; 
    DFmax = p.Results.dfmax;
    ylags = p.Results.ARlags;
    group_arlags = p.Results.GroupARlags;
    group_ctrl = p.Results.Groupctrl;
    gindex = p.Results.gindex;
    %nfolds = p.Results.nfolds; %If passed by user, it's ignored
    %foldid = p.Results.foldid; %If passed by user, it's ignored
    %nlambda = p.Results.nlambda;
    %lambda_factor = p.Results.lambda_factor;
    %lambda = p.Results.lambda;
    pf = p.Results.pf;
    pmax = p.Results.pmax;
    intercept = p.Results.intercept;
    eps = p.Results.eps;
    maxit = p.Results.maxit;
    peps = p.Results.peps;
    fe = p.Results.fe;
    Nfe = p.Results.N; 

    if isempty(oo) || isempty(fieldnames(oo)) %if the user supplies an empty struct OR empty array in 'Options'
        oo = struct('Method',{'CV-Kfold'}, 'Folds',10, 'OptimMeth',{'MinErr'});
    end

    valid_IC = {'BIC','AIC','AICc'};
    valid_CV = {'CV-Kfold'};
    validOptimMeth = [valid_IC, valid_CV];
    if ~any(strcmpi(oo.Method, validOptimMeth))
        error(['HPOptimizeOptions.Method can only be one of the following: ', char(join(validOptimMeth,", "))] )
    end
    
    %Default values for Name-Value pair 'HPOptimizeOptions' 
    if any(strcmpi(oo.Method,"CV-Kfold")) && ~isfield(oo,'Folds')
        oo.Folds = 10;
        warning('Name-Value pair ''HPOptimizeOptions'' did not contain field ''Folds''. ''Folds'' is set to 10.')
    end
    if any(strcmpi(oo.Method,valid_CV)) && ~isfield(oo,'OptimMeth')
        oo.OptimMeth = 'MinErr';
        warning('Name-Value pair ''HPOptimizeOptions'' did not contain field ''OptimMeth''. ''OptimMeth'' is set to ''MinErr''.')
    end
    if any(strcmpi(oo.Method,valid_CV))
        if strcmpi(oo.OptimMeth, 'MinErr')
            c=1;
        elseif strcmpi(oo.OptimMeth, 'MinErrPlus1SE')
            c=2;
        else
           error('Field ''OptimMeth'' of Name-Value pair ''HPOptimizeOptions'' contains invalid value. Choose between: MinErr, MinErrPlus1SE.') 
        end
    end

   if any(strcmpi(oo.Method,"CV-Kfold"))
       CVparam = oo.Folds;
   else
       CVparam = [];
   end
   
    alfa = sort(alfa(:),1,'descend'); %<- sort alfas descendingly
    
    [T,N] = size(X);
    Nd = size(D,2); %Number of variables in ctrl/D
    Ny = size(ylags,2); %Number of ARlags
    Xall = [X ylags D];
    Nall = size(Xall,2);
    
    %idxD = N+(1:Nd)';
    %idxY = N+Nd+(1:Ny)';
 
    if isempty(gindex)
        gindex = 1:N;
    end
    if ~(max(size(gindex))==N)
        warning('Name-Value pair ''gindex'' must be a vector with size(X,2) elements (i.e. excluding ''ARlags'' and ''ctrl'' variables). It has been set to default.')
        gindex = 1:N;
    end
    
    if isempty(D) %If user set group_ctrl=true, but did NOT provide D's 
        group_ctrl = false;
    end
    if isempty(ylags)
        group_arlags = false;
    end

    if group_ctrl
        gindex_D = ones(1,Nd);
    elseif ~group_ctrl 
        gindex_D = 1:Nd;
    end
    
    if group_arlags
        gindex_ar = ones(1,Ny);
    elseif ~group_arlags 
        gindex_ar = 1:Ny;
    end
    
    gindex = [gindex, max(gindex)+gindex_ar];
    gindex = [gindex, max(gindex)+gindex_D];
    
    if isempty(pf)
        pen = ones(Nall,1);
    else 
        pen = pf;
    end
    %pen([idxD idxY]) = 0; %Do NOT penalize Ylags & other control variables

   if any(strcmpi(oo.Method, valid_IC)) %optimal Lambda & Alfa via IC
       crit_a = nan(1, length(alfa)); 
       b_ica = nan(Nall+1, length(alfa));
       lambda_ica = nan(1, length(alfa)); 
       fit_ia = cell(1,length(alfa));
       for j = 1:length(alfa)
           fit = icsglfit(Xall,Y,oo.Method,'gamma',alfa(j),'gindex',gindex,'standardize',standardize,...
               'dfmax',DFmax,'pf',pen,'pmax',pmax,'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',fe,'N',Nfe);
           fit_ia{j} = fit;

           b_ica(:,j) = [fit.icsglfit.lam_min.b0; fit.icsglfit.lam_min.beta];
           crit_a(j) = min(fit.cvm);
           lambda_ica(j) = fit.lambda_min;
       end
       
       minICa = min(crit_a);
       idx = find(crit_a==minICa,1); 
       b = b_ica(:,idx);
       oV = crit_a(idx);
       oL = lambda_ica(idx); %Optimal lambda
       oA = alfa(idx);       %Optimal alfa
       fit = fit_ia{idx};
   
   elseif any(strcmpi(oo.Method, valid_CV)) %optimal Lambda & Alfa via CV
       %'c' -> 1=MinErr; 2=MinErrPlus1SE
       ca = c; 
       cl = c; 
       
       population = cat(2, repmat(1:CVparam, 1, floor(T/CVparam)), 1:mod(T,CVparam));
       cvfoldid = sort(population);

       pe_ia = nan(1,length(alfa));
       se_ia = nan(1,length(alfa));
       optCVidx_l = nan(2, length(alfa));
       optlambda_ia = nan(1,length(alfa));
       fit_ia = cell(1,length(alfa));
       for j = 1:length(alfa)
            fit = cvsglfit(Xall,Y,'gamma',alfa(j),'gindex',gindex,'foldid',cvfoldid,'parallel',parallel,'standardize',standardize,...
                'dfmax',DFmax,'pf',pen,'pmax',pmax,'intercept',intercept,'eps',eps,'maxit',maxit,'peps',peps,'fe',fe,'N',Nfe);
            fit_ia{j} = fit;
            lambdaseq = fit.sglfit.lambda;
           
            pe = fit.cvm; %CV-MSE for each LAMBDA (for the jth alfa)
            se  = fit.cvsd;

            minPE = min(pe);
            minIx = find(pe==minPE,1); 
            
            minplus1 = pe(minIx) + se(minIx);
            seIx = find((pe(1:minIx) <= minplus1),1,'first');

            optCVidx_l(:,j) = [minIx, seIx]'; 

            %CV-MSE of optimally selected Lambda (for given alfa)
            pe_ia(1,j) = pe(optCVidx_l(cl,j));
            se_ia(1,j) = se(optCVidx_l(cl,j));
            optlambda_ia(1,j) = lambdaseq(optCVidx_l(cl,j)); %Optimal Lambda
       end
       pe = pe_ia;
       se = se_ia;

       %1. AlfaMinMSE 
       minPE = min(pe);
       minIx = find(pe==minPE,1);

       %2. AlfaMinMSE1SE 
       minplus1 = pe(minIx) + se(minIx);
       seIx = find((pe(1:minIx) <= minplus1),1,'first');

       optCVidx_a = [minIx, seIx];
       opt.lambda = optlambda_ia(optCVidx_a(ca)); %Optimal Lambda
       opt.alfa = alfa(optCVidx_a(ca)); %Optimal Alfa
       
       %Extract betas from the 'sglfit' with the OPTIMAL alfa & lambda.
       %fit = sgl(Xall,Y,'gamma',opt.alfa,'gindex',gindex,'lambda',opt.lambda,...);
       fit = fit_ia{optCVidx_a(ca)}; %'sglfit' object corresponding to the optimal alfa
       lambdaseq = fit.sglfit.lambda;
       b_a = [fit.sglfit.b0'; fit.sglfit.beta];
       idx = find(lambdaseq==opt.lambda, 1,'first'); %the optimal lambda
       b = b_a(:,idx);
       oV = pe(optCVidx_a(ca));
       oL = lambdaseq(idx);
       oA = opt.alfa;
   end
   
   probj.b = b;
   probj.gindex = gindex;
   probj.OptimalModel.Value = oV; %The IC or CV-MSE of the selected model
   probj.OptimalModel.Lambda = oL; 
   probj.OptimalModel.Alfa = oA; 
   probj.OptimalModel.sgl = fit;
end