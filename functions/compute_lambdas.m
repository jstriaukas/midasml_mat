function lambdas = compute_lambdas(X, y, nlambda, lambda_min)
%--------------------------------------------------------------------------
% compute_lambdas: compute lambda sequence
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Compute lambda sequence based on data inputs and parameters.
%
% USAGE:
%    lambdas = compute_lambdas(X, y, nlambda, lambda_min)
%
% INPUT ARGUMENTS:
% x           Input matrix, of dimension nobs x nvars; each row is an
%             observation vector.
% y           Response variable. 
% nlambda     Number of lambdas values.
% lambda_min  Ratio of minimum to maximum lambda value.
%             
% OUTPUT ARGUMENTS:
% lambdas     The actual sequence of lambda values used.
%
% DETAILS:
%             The maximum lambda value is computed as 
%             lambda_max = |X'*y|_\infty. The minimum value
%             lam_min = lambda_max*lambda_min, where lambda_min
%             is the ratio of maximum and minimum lambda values.
%             If lambda_min is not supplied, lambda_min = 1e-2 for 
%             the case when nobs < nvars, else lambda_min = 1e-4.
%             lambdas = equidistanced points in the log space between 
%             lambda_max and lambda_min
%
% DATE: 2019-12-04
% AUTHOR: Jonas Striaukas
% LICENSE: GPL-2
%

[nobs, nvars] = size(X);
lam_max = max(abs(X'*y)/nobs);
if isempty(lambda_min)
    if nobs < nvars
        lambda_min = 1e-2;   
    else
        lambda_min = 1e-4;
    end
end
lam_min = lam_max*lambda_min;
log_lambdas = linspace(log(lam_max),log(lam_min), nlambda);
lambdas = exp(log_lambdas);



    