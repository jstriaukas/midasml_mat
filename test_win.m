clc;

addpath('functions');

load sglfitF_data

obj = sgl(x,y,'gamma',0.5);

save('win_solution', 'obj');


