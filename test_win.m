clc;

addpath('functions');

load sglfitF_data

obj = cvsglfit(x,y,'gamma',0.5);

save('win_solution', 'obj');


