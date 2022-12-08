%--------------------------------------------------------------------------
%
%   This file provides an example of how to apply sg-LASSO-MIDAS model.
%   
%   I use the famous FRED-MD dataset which is available at
%   https://research.stlouisfed.org/econ/mccracken/fred-databases/
%   I use the file 2021-10.csv downloaded on 2021/11/28 (and data cleaning
%   files) 
%
%--------------------------------------------------------------------------

clear; clc;

addpath '../functions';

% data

fred_codes = ["GDPC1";"CONSUMER";"RPI";"INDPRO";"IPMANSICS";"TCU";"PAYEMS";"UNRATE";"HOUST";"CPIAUCSL"];
[x, Meta] = fred_api(fred_codes, 'StartDate',datetime(1995,1,1));
data = x{:,:};

% let's take the first var in FRED MD and predict it using all other vars
% in the database using 12 monthly lags

y = data(:,1);
xtmp = data(:,2:end);

% data y and xtmp contains NaNs. Some are in the middle of the sample
% (perhaps missing observations) some are at the end (most likely ragged
% edge stuff). I impute zeros for simplicity. To deal with ragged edge
% problem, you can shift the data accordingly. 

y(isnan(y)) = 0;
xtmp(isnan(xtmp)) = 0;

% number of lags:
jmax = 12;
horizon = 1;
% example of data:
%mixed_freq_data(xtmp(:,1), dates, xtmp(:,1), dates, jmax, horizon, 0, dates(1), dates(end), true);


% choose degree fpr Legendre polynomials:
degree = 3; % 3 is typically good enough to be able to capture nonlinearities of lag/MIDAS weight function.

g = linspace(0, 1, jmax)';
w = lb(g, degree, 0, 1)/jmax;% weights of Legendre poly.

warning('off')

xtmp = xtmp(397:end,:);
y = y(397:end);

[T,p] = size(xtmp);
x = zeros(T - jmax - horizon + 2, p * (degree + 1));
gindex = zeros(p * (degree + 1),1);

for i = 1:p
   xtmpn = xtmp(:,i);
   out = mixed_freq_data(xtmpn,dates,xtmpn,dates,jmax,horizon,0,dates(1),dates(end),false);
   xw = out.EstX*w;
   idx = (i * (degree + 1) - degree) : i * (degree + 1);
   x(:,idx) = xw;
   gindex(idx) = repmat(i, degree + 1, 1);
end

% drop data from the start of the y vector to match x (due to lags etc.):
y = y((jmax + horizon - 1):end);


% run cv sg-LASSO:
% takes some time...
fit = cvsglfit(x, y, 'gamma', 0.5, 'gindex', gindex, 'standardize', true);
fit = cvsglfit(x, y, 'gamma', 0.4, 'gindex', gindex, 'standardize', true);

% beta's:
fit.cvsglfit.lam_min.beta
% intercept:
fit.cvsglfit.lam_min.b0


% compute CV solution for both gamma and lambda:
% run loop over gamma grid (takes a bit more time....).
% gamma = linspace(1,0,21);
% cvms = zeros(21, 1);
% for gi = 1:21
%     fit = cvsglfit(x, y, 'gamma', gamma(gi), 'gindex', gindex);
%     cvms(gi) = min(fit.cvm);
% end
% idx = find(min(cvms)==cvms,1,'first');
% % optimal gamma:
% gamma(idx)
% fit = cvsglfit(x, y, 'gamma', gamma(idx), 'gindex', gindex);
