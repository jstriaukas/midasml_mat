function [pred,beta] = icsglfitpredict(obj, newx, varargin)       

%SYNTAX: 
%[yhat,bz] = icsglfitpredict(obj, newX, fe)

%OPTIONAL INPUTS:
%- 'fe' set to 'true' for Panel data (fixed effects), set to 'false' otherwise.
%       Default: fe=false

%Examples:
%[yhat,bz] = icsglfitpredict(obj,newX)
%[yhat,bz] = icsglfitpredict(obj,newX,false) %Equivalent to above
%[yhat,bz] = icsglfitpredict(obj,newX,true)  %Prediction for Panel/Fixed Effects

p = inputParser;
addRequired(p,'obj',@(z) isstruct(z));
addRequired(p,'newx',@isnumeric);
addOptional(p,'fe',false,@(z) isnumeric(z) || islogical(z));

parse(p, obj, newx,varargin{:});
obj = p.Results.obj;
newx = p.Results.newx;
fe = p.Results.fe;

[pred,beta] = cvsglfitpredict(obj,newx,'lambda_min',fe);

end
