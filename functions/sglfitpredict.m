function pred = sglfitpredict(obj, newx, fe)
[n,~] = size(newx);
if fe
    fes = [];
    for i = 1:obj.N
        fes = [fes; repmat(obj.a0(i,:),n/obj.N,1)];
    end
    pred = newx*obj.beta + fes;
    else
    pred = newx*obj.beta + repmat(obj.b0', n, 1);
end
end
