clear all
global xx
format long
xx = 0;
tic

%     [factor,fval] = ga(@calerror,3,[],[],[],[],[1.2042e-9 0.0439 0.002], [1.4542e-9 0.0728 0.005 ])
    PSO(@calerror,[],struct('Xdim',3,'Xmin',[1.2042e-9 0.0439 0.002]','Xmax',[1.4542e-9 0.0728 0.005]'))

toc



% factor =
% 
%    0.000000001220795   0.067847036036459   0.002028615779498
% 
% 
% ans =
% 
%   -0.347117868136479
