function [Psi] = lb(X, degree, a, b)
% lb Legendre polynomials shifted to [a,b]
%   For a given set of points in X, computes the orthonormal Legendre 
%   polynomials basis of L2[a,b] for a given degree
%   Returns a dim(X) x (degree+1) matrix
%
%   Required arguments:
%   X         = n x 1 vector of evaluation points
%   degree    = 1 x 1 value of the degree of polynomials
%   a and b   = 1 x 1 value each for [a, b] interval defining L2[a,b]. By
%   defult a  = b = 0
%
%   References:
%   H. J. Weber, G. B. Arfken, Essential Mathematical Methods for Physicists,
%   Elsevier Academic Press, San Diego, USA, 2004.

[n, ~] = size(X);

if nargin < 3
    a = 0;
    b = 1;
end

P = ones(n, degree+2); 
P(:, 2) = 2*X/(b-a) - (b+a)/(b-a);
Psi = ones(n, degree+1) / sqrt(b-a);

for i = 1:degree
    P(:, i+2) = (2*i+1)/(i+1) * P(:, 2) .* P(:, i+1) - i/(i+1) * P(:, i);
    Psi(:, i+1) = sqrt((2*i + 1) / (b-a)) * P(:, i+1);
end

end