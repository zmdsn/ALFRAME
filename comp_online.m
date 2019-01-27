clear all;
clc;  
format short g;
addroad;
totaltime = cputime;   
rng('default');

% **************************==- func_config -==****************************
func_name = {};
func_config = [];
% cec15 
func_name = [func_name {'cec15_f1','cec15_f2','cec15_f3','cec15_f4','cec15_f5','cec15_f6','cec15_f7','cec15_f8','cec15_f9','cec15_f10'}];
func_name = [func_name {'cec15_f11','cec15_f12','cec15_f13','cec15_f14','cec15_f15'}];

% cec15nich
% func_name = [func_name {'cec15ni_f1','cec15ni_f2','cec15ni_f3','cec15ni_f4','cec15ni_f5','cec15ni_f6','cec15ni_f7','cec15ni_f8'}];
% func_name = [func_name {'cec15ni_f9','cec15ni_f10','cec15ni_f11','cec15ni_f12','cec15ni_f13','cec15ni_f14','cec15ni_f15'}];

% cec15exp
% func_name = [func_name {'cec15exp_f1','cec15exp_f2','cec15exp_f3','cec15exp_f4','cec15exp_f5','cec15exp_f6','cec15exp_f7','cec15exp_f8'}];
% func_name = [func_name {'cec15exp_f9','cec15exp_f10','cec15exp_f11','cec15exp_f12','cec15exp_f13','cec15exp_f14','cec15exp_f15'}];

% old functions
% func_name = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'};
% func_name = {'CHf1','CHf2'};

% cec14
% func_name = [func_name {'cec14_f1','cec14_f2','cec14_f3','cec14_f4','cec14_f5','cec14_f6','cec14_f7','cec14_f8','cec14_f9','cec14_f10'}];
% func_name = [func_name {'cec14_f11','cec14_f12','cec14_f13','cec14_f14','cec14_f15','cec14_f16','cec14_f17','cec14_f18','cec14_f19','cec14_f20'}];
% func_name = [func_name {'cec14_f21','cec14_f22','cec14_f23','cec14_f24','cec14_f25','cec14_f26','cec14_f27','cec14_f28','cec14_f29','cec14_f30'}];

% func_name = [func_name {'f3'}];
% func_name = {'f1','f2','f3','f6','Cf1','Cf2','Cf3','Cf6'};

nFunc_num = length(func_name);
func_config.Xdim = 30;

% *************************==- algo_config -==*****************************
algo_name = { 'PSO','OLPSO','CLPSO','SPSO2011'};
% algo_name = { 'HS','NGHS','SGHS','IGHS','IHS'};
% algo_name = { 'bRGA'};
% algo_name = {'DErand1','SaDE','CoDE','GAD'};

% algo_name = { 'PSO','CLPSO','OLPSO','CEPSO'};
% algo_name = { 'PSO','CLPSO','OLPSO','DErand1','SE1','SE2','SE','bRGA'};
% algo_name = { 'SE2','SE','bRGA'};,'SE2','SE',,'SaDE','OLPSO'

% algo_name = { 'PSO','CLPSO','DErand1','bRGA','jDE','SE1'};

% algo_name = { 'SE','SE3','SE31','SE12','SE1'};,'SPSO2011''SE','CoDE',
% algo_name = { 'ODE','OCDE','ODDE'};
% algo_name = { 'LTES','LTXES','CGA','PSO','DErand1', 'bRGA'};
% algo_name = { 'PSO','CLPSO','OLPSO','DErand1', 'bRGA', 'RGA','HAF1','HAF2','HAF3'};
% algo_name = {  'PSO','CLPSO','OLPSO','DErand1', 'RGA','HAF1','HAF2','HAF3'};
% algo_name = { 'PSO', 'HAF3'};

algo_config.PopuSize = 100;
algo_config.maxFES = func_config.Xdim *10000; %2* 

nAlgorithm = length(algo_name);

% *************************==- Frame_config -==*****************************
Run = 5;             
result_folder = 'TesT/'; 


finalRes = zeros(nFunc_num*nAlgorithm,5);
Result_val = zeros(5);
Algorithm = cell(nAlgorithm,Run);
check_frame; 
for func_num = 1:nFunc_num
    fxopt = zeros(1,nAlgorithm); 
    Fx = zeros(Run,nAlgorithm); 
    for j=1:Run
      display(['SRS_TEST F' num2str(func_num) ', RUN' num2str(j) ' begins ==>>']);
      for  k = 1:nAlgorithm
        [fxopt(k),best_solution, Algorithm{k}{j}, run_info]= feval([algo_name{k}],[func_name{func_num}],algo_config,func_config);
%         disp( ['  -==********==- ' algo_name{k} ' -==*********==-'] );
%         disp( ['      Best Value  =  ' num2str(fxopt(k),'%10.5e') ]);
%         disp( ['      Time Used   =  ' num2str(run_info) 's']);
      end
      Fx(j,:) =  fxopt;
    end
    
    % **************************==- analysis -==***************************
    % ********==- draw online figure -==********
    draw_data = cell(nAlgorithm);
    for  k = 1:nAlgorithm
      draw_data{k} = [];
      
      
      draw_data{k}(:,1) = Algorithm{k}{1}(:,1);
      draw_data{k}(:,2) = Algorithm{k}{1}(:,2);
      for j = 2 : Run
        l1 = length(draw_data{k}(:,1));
        l2 = length(Algorithm{k}{j}(:,2));
        if l1 == l2
            draw_data{k}(:,2) = draw_data{k}(:,2) + Algorithm{k}{j}(:,2);
        else
            if l1 > l2
                draw_data{k}(:,2) = draw_data{k}(:,2) + [Algorithm{k}{j}(:,2);ones(l1 - l2,1)*Algorithm{k}{j}(end,2) ];
            else 
                draw_data{k}(:,2) = draw_data{k}(:,2) + Algorithm{k}{j}(1:l1,2);
            end
        end
      end
      leng_str{k} = ['\fontsize{15}' algo_name{k}];
    end
    close all;
    draw_lines(draw_data , nAlgorithm, Run , func_num ,leng_str,1,result_folder);
    
    
    % ********==- Wilcoxon -==********
    Result_val = [min(Fx)',mean(Fx)',median(Fx)',std(Fx)'];
    for j=1:nAlgorithm     % ttest
        [P,H,STATS] = ranksum(Fx(:,j), Fx(:,nAlgorithm));
        finalRes((func_num-1)*nAlgorithm+j,:) = [Result_val(j,:) sign(finalRes(func_num*nAlgorithm,2) - finalRes((func_num-1)*nAlgorithm+j,2))*H];
    end
    
    RSS = [];
    for  j = 1:nAlgorithm
        k = (func_num-1)*nAlgorithm+j;
        RSS = [RSS finalRes(k,2) finalRes(k,4) finalRes(k,5)];
    end   
    Sr(func_num,:) = RSS;
    
    % ********==- draw boxfig -==********
    close all;
    boxplot(Fx);
    set(gca,'xticklabel',algo_name);
    saveas( gcf , [result_folder 'boxfig/' func_name{func_num}  '_box.fig']);
    save( [result_folder 'data/' func_name{func_num}  '_result.mat'] , 'Fx','draw_data');
end  %endfor func_name



% ********==- data carding -==********
Transpose = [];
for  j = 1:nFunc_num
    k1 = (j-1)*nAlgorithm + 1;
    k2 = j*nAlgorithm;
    Transpose = [Transpose ; finalRes(k1:k2,:)'];
end   

cal_Sr = zeros(4,nAlgorithm);
for j=1:nAlgorithm     % ttest
    cal_Sr(1,j) = sum(Sr(:,3*j) == 1);
    cal_Sr(2,j) = sum(Sr(:,3*j) == -1);
    cal_Sr(3,j) = sum(Sr(:,3*j) == 0);
    cal_Sr(4,j) = cal_Sr(2,j) - cal_Sr(1,j);
end

save( [result_folder 'result.mat'] , 'finalRes','Sr','Transpose','cal_Sr');


format short;
display(['Total time ',num2str(cputime - totaltime)]);








