function [pred,beta] = cvsglfitpredict(obj, newx, varargin)       

%SYNTAX: 
%[yhat,bz] = cvsglfitpredict(obj, newX, s, fe)

%OPTIONAL INPUTS:
%- 's' is the choice of the optimal lambda: 'lambda_min' OR 'lambda_1se'
%      Default: s=lambda_min
%- 'fe' set to 'true' for Panel data (fixed effects), set to 'false' otherwise.
%       Default: fe=false

%Examples:
%[yhat,bz] = cvsglfitpredict(obj,newX) 
%[yhat,bz] = cvsglfitpredict(obj,newX,'lambda_min') %Equivalent to above
%[yhat,bz] = cvsglfitpredict(obj,newX,'lambda_min',false) %Equivalent to above
%[yhat,bz] = cvsglfitpredict(obj,newX,'lambda_1se')

p = inputParser;
addRequired(p,'obj',@(z) isstruct(z));
addRequired(p,'newx',@isnumeric);
addOptional(p,'s','lambda_min',@(z) isstring(z) || ischar(z) || iscell(z));
addOptional(p,'fe',false,@(z) isnumeric(z) || islogical(z)); %Optional input. Default: fe=false
parse(p, obj, newx,varargin{:});
obj = p.Results.obj;
newx = p.Results.newx;
s = p.Results.s;
fe = p.Results.fe;

if strcmp(s,'lambda_min')
    opt_lam = obj.lambda_min;
elseif strcmp(s,'lambda_1se')
    opt_lam = obj.lambda_1se;
end 

idx = find(obj.sglfit.lambda==opt_lam);

newobj.b0 = obj.sglfit.b0(idx);
newobj.beta = obj.sglfit.beta(:,idx);
if fe 
    newobj.N = obj.sglfit.N; 
    newobj.a0 = obj.sglfit.a0(:,idx);
end
[pred,beta] = sglfitpredict(newobj, newx, fe);

end




