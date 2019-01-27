% Main function
clear all;
global  initial_flag
initial_flag = 0;
format long e;
D = 30;%10
NP =1000;
Xmin1=0*ones(1,D);
Xmax1=+10*ones(1,D);
Xmin2=-5.12*ones(1,D);
Xmax2=+5.12*ones(1,D);
Xmin3=-1000*ones(1,D);
Xmax3=+1000*ones(1,D);
Xmin4=-50*ones(1,D);
Xmax4=+50*ones(1,D);
Xmin5=-600*ones(1,D);
Xmax5=+600*ones(1,D);
Xmin6=-600*ones(1,D);
Xmax6=+600*ones(1,D);
Xmin7=-140*ones(1,D);
Xmax7=+140*ones(1,D);
Xmin8=-140*ones(1,D);
Xmax8=+140*ones(1,D);
Xmin9=-500*ones(1,D);
Xmax9=+500*ones(1,D);
Xmin10=-500*ones(1,D);
Xmax10=+500*ones(1,D);
Xmin11=-100*ones(1,D);
Xmax11=+100*ones(1,D);
Xmin12=-1000*ones(1,D);
Xmax12=+1000*ones(1,D);
Xmin13=-500*ones(1,D);
Xmax13=+500*ones(1,D);
Xmin14=-1000*ones(1,D);
Xmax14=+1000*ones(1,D);
Xmin15=-1000*ones(1,D);
Xmax15=+1000*ones(1,D);
Xmin16=-10*ones(1,D);
Xmax16=+10*ones(1,D);
Xmin17=-10*ones(1,D);
Xmax17=+10*ones(1,D);
Xmin18=-50*ones(1,D);
Xmax18=+50*ones(1,D);
Max_FES=D*20000;
Max_Gen=D*100;
DELTA=10^(-4);

gn=[2 2 0 0 0 0 0 0 0 0 1 1 3 3 3 2 2 1];
hn=[0 1 1 4 2 2 1 1 1 1 0 0 0 0 0 2 1 1];

func_num=18; % Function Number
    
eval(['Xmin=Xmin' int2str(func_num) ';']);
eval(['Xmax=Xmax' int2str(func_num) ';']); 
    
pop=repmat(Xmin,NP,1)+repmat((Xmax-Xmin),NP,1).*rand(NP,D);
size(pop)
[val(:,1), g, h] = fcnsuite_func(pop,func_num);


   