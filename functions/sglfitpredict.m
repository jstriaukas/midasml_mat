function [pred,beta] = sglfitpredict(obj, newx, fe)
if nargin < 3
    fe = false;
end

[n,~] = size(newx);
if fe
    fes = [];
    b0 = obj.a0;
    for i = 1:obj.N 
        fes = [fes; repmat(b0(i,:), n/obj.N,1)];
    end
    pred = newx*obj.beta + fes;
    
else
    b0 = obj.b0';
    pred = newx*obj.beta + repmat(b0, n, 1); %size(newX,1)-by-L
end

%obj.beta is size(x,2)-by-L
%b0 is: 
%- 1xL for fe=false => beta is size(x,2)+1-by-L
%- NxL for fe=true => beta is [size(x,2)+N]-by-L
beta = [b0; obj.beta]; 

end
