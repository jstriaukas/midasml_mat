%--------------------------------------------------------------------------
%
%   This file computes the data for the example of how to apply 
%   sg-LASSO-MIDAS model.
%   
%   I use the famous FRED-MD dataset which is available at
%   https://research.stlouisfed.org/econ/mccracken/fred-databases/
%   I use the file 2021-10.csv downloaded on 2021/11/28 (and data cleaning
%   files) 
%   Data cleaning files were downloaded from the same page. 
%   I used fredfactors.m upto removing outliers line 111 in the original
%   code. Be aware that such procedures may significantly impact your
%   results (see comment below). 
%
%--------------------------------------------------------------------------

clear
close all
clc

% =========================================================================
% DESCRIPTION 
% This script loads in a FRED-MD dataset, processes the dataset, and then
% estimates factors.
%
% -------------------------------------------------------------------------
% BREAKDOWN OF THE SCRIPT
% 
% Part 1: Load and label FRED-MD data.
%
% Part 2: Process data -- transform each series to be stationary and remove
%         outliers.
%
% Part 3: Estimate factors and compute R-squared and marginal R-squared. 
%
% -------------------------------------------------------------------------
% AUXILIARY FUNCTIONS
% List of auxiliary functions to be saved in same folder as this script.
%
%   prepare_missing() - transforms series based on given transformation
%       numbers
%
%   remove_outliers() - removes outliers
%
%   factors_em() - estimates factors
%
%   mrsq() - computes R-squared and marginal R-squared from factor 
%       estimates and factor loadings
%
% -------------------------------------------------------------------------
% NOTES
% Authors: Michael W. McCracken and Serena Ng
% Date: 9/5/2017
% Version: MATLAB 2014a
% Required Toolboxes: None
%
% -------------------------------------------------------------------------
% PARAMETERS TO BE CHANGED

% File name of desired FRED-MD vintage
csv_in='2021-10.csv';

% Type of transformation performed on each series before factors are
% estimated
%   0 --> no transformation
%   1 --> demean only
%   2 --> demean and standardize
%   3 --> recursively demean and then standardize
DEMEAN=2;

% Information criterion used to select number of factors; for more details,
% see auxiliary function factors_em()
%   1 --> information criterion PC_p1
%   2 --> information criterion PC_p2
%   3 --> information criterion PC_p3
jj=2;

% Maximum number of factors to be estimated; if set to 99, the number of
% factors selected is forced to equal 8
kmax=8;

% =========================================================================
% PART 1: LOAD AND LABEL DATA

% Load data from CSV file
dum=importdata(csv_in,',');

% Variable names
series=dum.textdata(1,2:end);

% Transformation numbers
tcode=dum.data(1,:);

% Raw data
rawdata=dum.data(2:end,:);

% JS: modified original code to get proper matlab dates....
dates = datenum(dum.textdata(3:end,1));

% T = number of months in sample
T=size(dates,1);
rawdata=rawdata(1:T,:);

% =========================================================================
% PART 2: PROCESS DATA

% Transform raw data to be stationary using auxiliary function
% prepare_missing()
yt=prepare_missing(rawdata,tcode);

% Reduce sample to usable dates: remove first two months because some
% series have been first differenced
yt=yt(3:T,:);
dates=dates(3:T,:);

%--------------------------------------------------------------------------
% JS: CAUTION!!! The function is taken as it is from 
% https://research.stlouisfed.org/econ/mccracken/fred-databases/
%  1. Results are inevitably sensitive to the method of removing outliers 
%  2. Empirical risk minimization may be affected by the presence of
%     outliers in BOTH response and covariates
%--------------------------------------------------------------------------

% Remove outliers using auxiliary function remove_outliers(); see function
% or readme.txt for definition of outliers
%   data = matrix of transformed series with outliers removed
%   n = number of outliers removed from each series
[data,n]=remove_outliers(yt);


save('fred_example_data.mat', 'data', 'dates');