clc;

addpath('functions');

load sglfitF_data

obj1 = sgl(x,y, 'gamma', 0.5);

clearvars -except obj1

load win_solution.mat

if ~all(obj1.b0==obj.b0)
    error('b0 not the same');
end

if ~all(obj1.beta==obj.beta)
    error('b0 not the same');
end

if ~all(obj1.npass==obj.npass)
    error('b0 not the same');
end

if ~all(obj1.lambda==obj.lambda)
    error('b0 not the same');
end



