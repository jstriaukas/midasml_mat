function [pred,beta] = cvsglfitpredict(obj, newx, s)       
    
if nargin < 3
    s = 'lambda_min';
end

if strcmp(s,'lambda_min')
    opt_lam = obj.lambda_min;
elseif strcmp(s,'lambda_1se')
    opt_lam = obj.lambda_1se;
end 

idx = find(obj.sglfit.lambda==opt_lam);
b0 = obj.sglfit.b0(idx);
beta =obj.sglfit.beta(:,idx);
pred = newx*beta+b0;
beta = [b0;beta];

end




