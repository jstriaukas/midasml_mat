load 'fred_example/fred_example_data.mat';
addpath 'functions';
y = data(:,1);
xtmp = data(:,2:end);
y(isnan(y)) = 0;
xtmp(isnan(xtmp)) = 0;
jmax = 12; % number of high-frequency lags
horizon = 1;
degree = 3;
g = linspace(0,1, jmax)';
a=0; b=1;
w = lb(g, degree, a,b)/jmax;% WEIGHTS of Legendre polynomial 
[T, p] = size(xtmp);

x = zeros(T-jmax-horizon+2, p*(degree+1));
gindex = zeros(p*(degree+1), 1);
warning('off') 
for i = 1:p %loop over 126 monthly variables
   xtmpn = xtmp(:,i);
   out = mixed_freq_data(xtmpn,dates,xtmpn,dates,jmax,horizon,0,dates(1),dates(end),false);
   xw = out.EstX*w;
   idx = (i*(degree+1)-degree):i*(degree+1);
   x(:, idx) = xw;
   gindex(idx) = repmat(i, degree+1, 1);
end

y = y((jmax+horizon-1):end);
clearvars -except x y gindex newX

%Example 1: This works perfectly!
o = sgl(x,y,'gamma',0.5,'lambda',[0.1, 0.0005, 1e-04],'gindex',gindex,'standardize',true);

%Example 2: This produces INFINITE LOOP!
o = sgl(x,y,'gamma',0.5,'lambda',[],'gindex',gindex,'standardize',true); 