function [Psi] = gb(X, alpha, degree, a, b)
% gb Gegenbauer polynomials shifted to [a,b]
[n, ~] = size(X);
if nargin < 3
    a = 0;
    b = 1;
end

P = ones(n, degree+2); 
P(:, 2) = 2*X/(b-a) - (b+a)/(b-a);
Psi = ones(n, degree+1) / sqrt(b-a);

for i = 1:degree
    d = (2 * i + 2 * alpha + 1) * (2 * i + 2 * alpha + 2)/(2 * (i + 1) * (i + 2 * alpha + 1));
    c = (alpha + i)^2 * (2 * i + 2 * alpha + 2)/((i + 1) * (i + 2 * alpha + 1) * (2 * i + 2 * alpha));
    P(:, i+2) = d .* P(:, 2) .* P(:, i + 1) - c  .* P(:, i);
    Psi(:, i+1) = sqrt((2 * i + 1)/(b - a)) * P(:, i + 1);
end

end
