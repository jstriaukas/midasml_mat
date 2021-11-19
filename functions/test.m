load sglfitF_data

obj = sgl(x,y);

obj = cvsglfit(x,y);

obj = icsglfit(x,y);

obj = sgl(x,y, 'fe', true, 'N', 10);
