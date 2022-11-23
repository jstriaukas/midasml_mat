function ulam = computelambda(nlambda, flmin, nobs, x, y, gamma, gindex, ngroups, pf, standardize)
    maxlam = maxlambda(nobs, x, y, gamma, gindex, ngroups, pf, standardize);
    if isnan(maxlam)
        error('data in x or y has missing entries/NaNs. program cannot proceed');
    end
    ulam = double(zeros(nlambda,1));
    ulam(1) = maxlam;
    for j = 2:nlambda
                tmp = log(maxlam) + (log(maxlam*flmin) - log(maxlam)) * (j - 1) / (nlambda - 1);
                ulam(j) = exp(tmp);
    end
end


function maxlam = maxlambda(nobs, X, y, gamma, gindex, ngroups, pf, standardize)
    r = y;
    x = X;
    if standardize
        r = zscore(r);
        x = zscore(x);
    end 
    xy = x'*r;  
    xy = xy/double(nobs);
    wmaxg = zeros(ngroups,1);
    
    if (gamma == 1.0) 
        maxlam = max(abs(xy));
    else 
       for k = 1:ngroups
           gend = gindex(k);
           if k == 1
               gstart = 1;
           else
               gstart = gindex(k-1) + 1;
           end
           gw = 0;
           for gj = gstart:gend
               gw = gw + pf(gj);
           end
           gw = sqrt(gw);
           if gw == 0.0 
                wmaxg(k) = 0.0;
           else
               if gamma == 0.0 
                   rb = sqrt(xy(gstart:gend)'*xy(gstart:gend));
                   wmaxg(k) = rb/gw;
               else
                   lb = 0;
                   rb = max(abs(xy(gstart:gend)))/gamma;
                   rb = solvewmaxg(gstart, gend, gamma, lb, rb, gw, pf, xy);
                   wmaxg(k) = rb;
               end
           end
       end
       maxlam = max(wmaxg);
    end
end

function out = solvewmaxg(gstart, gend, gamma, lb, rb, gw, pf, xy)
    stopflag = false;
    tol = 1e-13;
    while stopflag == false
        mp = 0.5 * (lb + rb);
        fl = 0.0;
        fm = 0.0;
        fr = 0.0;
        for indexi =  gstart:gend
            tmpl = abs(xy(indexi)) - gamma * lb * pf(indexi);
            tmpm = abs(xy(indexi)) - gamma * mp * pf(indexi);
            tmpr = abs(xy(indexi)) - gamma * rb * pf(indexi);
            if tmpl > 0.0
                fl = fl + tmpl * tmpl;
            end
            if tmpm > 0.0
                fm = fm + tmpm * tmpm;
            end
            if tmpr > 0.0
                fr = fr + tmpr * tmpr;
            end
        end
        fl = fl - (1.0 - gamma) * (1.0 - gamma) * lb * lb * gw * gw;
        fm = fm - (1.0 - gamma) * (1.0 - gamma) * mp * mp * gw * gw;
        fr = fr - (1.0 - gamma) * (1.0 - gamma) * rb * rb * gw * gw;
        if fl * fm < 0.0
            if abs(lb - mp) > tol
                rb = mp;
            else
                stopflag = true;
            end
        else
            if fm * fr < 0.0
                if abs(mp - rb) > tol
                    lb = mp;
                else
                    stopflag = true;
                end
            else
                stopflag = true;
            end
        end
    end
    out = mp;
end