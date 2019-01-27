function [f,g,h]=fcnsuite_func(x,func)
x=x';
if func==1
    [f, g, h] = mlbsuite(x, 2, 0, 'C01');
elseif func==2
    [f, g, h] = mlbsuite(x, 2, 1, 'C02');
elseif func==3
    [f, g, h] = mlbsuite(x, 0, 1, 'C03');
elseif func==4
    [f, g, h] = mlbsuite(x, 0, 4, 'C04');
elseif func==5
    [f, g, h] = mlbsuite(x, 0, 2, 'C05');
elseif func==6
    [f, g, h] = mlbsuite(x, 0, 2, 'C06');
elseif func==7
    [f, g, h] = mlbsuite(x, 1, 0, 'C07');
elseif func==8
    [f, g, h] = mlbsuite(x, 1, 0, 'C08');
elseif func==9
    [f, g, h] = mlbsuite(x, 0, 1, 'C09');
elseif func==10
    [f, g, h] = mlbsuite(x, 0, 1, 'C10');
elseif func==11
    [f, g, h] = mlbsuite(x, 0, 1, 'C11');
elseif func==12
    [f, g, h] = mlbsuite(x, 1, 1, 'C12');
elseif func==13
    [f, g, h] = mlbsuite(x, 3, 0, 'C13');
elseif func==14
    [f, g, h] = mlbsuite(x, 3, 0, 'C14');
elseif func==15
    [f, g, h] = mlbsuite(x, 3, 0, 'C15');
elseif func==16
    [f, g, h] = mlbsuite(x, 2, 2, 'C16');
elseif func==17
    [f, g, h] = mlbsuite(x, 2, 1, 'C17');
elseif func==18
    [f, g, h] = mlbsuite(x, 1, 1, 'C18');
end