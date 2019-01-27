clear all;
clc;  
format short g;
addroad;
totaltime = cputime;   
rng('default')

% **************************==- func_config -==****************************
func_name = {};
func_config = [];
% func_name = [func_name {'cec15_f2'}];

% cec15  
func_name = [func_name {'cec15_f1','cec15_f2','cec15_f3','cec15_f4','cec15_f5','cec15_f6','cec15_f7','cec15_f8','cec15_f9','cec15_f10'}];
func_name = [func_name {'cec15_f11','cec15_f12','cec15_f13','cec15_f14','cec15_f15'}];
% % func_name = [func_name {'cec15_f7','cec15_f8','cec15_f9','cec15_f10'}];
% % func_name = [func_name {'cec15_f13','cec15_f14','cec15_f15'}];
% % func_name = [func_name {'cec15_f5','cec15_f6','cec15_f15'}];
% % func_name = [func_name {'cec15_f1','cec15_f6','cec15_f8','cec15_f10','cec15_f14'}];
% % func_name = [func_name {'cec15_f14'}];


% % cec15nich
% func_name = [func_name {'cec15ni_f1','cec15ni_f2','cec15ni_f3','cec15ni_f4','cec15ni_f5','cec15ni_f6','cec15ni_f7','cec15ni_f8'}];
% func_name = [func_name {'cec15ni_f9','cec15ni_f10','cec15ni_f11','cec15ni_f12','cec15ni_f13','cec15ni_f14','cec15ni_f15'}];

% % cec15exp
% func_name = [func_name {'cec15exp_f1','cec15exp_f2','cec15exp_f3','cec15exp_f4','cec15exp_f5','cec15exp_f6','cec15exp_f7','cec15exp_f8'}];
% func_name = [func_name {'cec15exp_f9','cec15exp_f10','cec15exp_f11','cec15exp_f12','cec15exp_f13','cec15exp_f14','cec15exp_f15'}];

% old functions
% func_name = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'};
% func_name = {'CHf1','CHf2'};

% cec14
% func_name = [func_name {'cec14_f1','cec14_f2','cec14_f3','cec14_f4','cec14_f5','cec14_f6','cec14_f7','cec14_f8','cec14_f9','cec14_f10'}];
% func_name = [func_name {'cec14_f11','cec14_f12','cec14_f13','cec14_f14','cec14_f15','cec14_f16','cec14_f17','cec14_f18','cec14_f19','cec14_f20'}];
% func_name = [func_name {'cec14_f21','cec14_f22','cec14_f23','cec14_f24','cec14_f25','cec14_f26','cec14_f27','cec14_f28','cec14_f29','cec14_f30'}];
nFunc_num = length(func_name);
func_config.Xdim = 30; % 10 30 50 100


% *************************==- algo_config -==*****************************
algo_name = [];
algo_name = { 'PSO','OLPSO','CLPSO','SPSO2011','SPSO2007'};
% algo_name = [algo_name { 'HS','NGHS','SGHS','IHS'}];
algo_name = [algo_name { 'bRGA'}];
algo_name = [algo_name {'DErand1','DErand2','DEbest1','DEbest2','DEctob1','SaDE','CoDE','JADE','jDE'}];
algo_name = [algo_name { 'CMAES'}];

% algo_name = {  'HAF11','HAF21', 'HAF31', 'HAF12','HAF22', 'HAF32', 'HAF13','HAF23', 'HAF33', 'HAF14','HAF24', 'HAF34', 'HAF15','HAF25', 'HAF35'};

% algo_name = [algo_name { 'CMAES','GAD'}];
% algo_name = {  'PSO','DErand1', 'RGA','HAF1','HAF2','HAF3'};
% algo_name = {  'HAF21','HAF22', 'HAF23','HAF24','HAF25', 'HAF26','HAF27'};
% algo_name = [algo_name {  'HAF31','HAF32', 'HAF33','HAF34','HAF35', 'HAF36','HAF37'}];

% algo_name = { 'PSO','CLPSO','OLPSO','CEPSO'};,'IGHS','GAD','APSO'
% algo_name = { 'PSO','CLPSO','OLPSO','DErand1','SE1','SE2','SE','bRGA'};
% algo_name = { 'SE2','SE','bRGA'};,'SE2','SE',,'SaDE','OLPSO','SE1'

% algo_name = { 'PSO','CLPSO','DErand1','bRGA','jDE'};
% algo_name = [algo_name { 'SE'}];


algo_config.PopuSize = 100;
algo_config.maxFES = func_config.Xdim * 10000;

nAlgorithm = length(algo_name);

% *************************==- Frame_config -==*****************************

Run = 30;                   % runtimes
totle_RUN = 60;
result_folder = 'Archives/';    % the folder to save the results


% **************************==- main loop -==******************************

finalRes = zeros(nFunc_num*nAlgorithm,5);
Result_val = zeros(5);
Algorithm = cell(nAlgorithm,Run);
check_frame; % check_frame
for func_num = 1:nFunc_num
    % **************************==- alrecord -==***************************
    fxopt = zeros(1,nAlgorithm); 
    Fx = zeros(Run,nAlgorithm);  
    display(['SRS_TEST F' num2str(func_num) ',  begins ==>>']);
    for  k = 1:nAlgorithm
        dir_folds = [result_folder  algo_name{k} '/'];
        if ~exist( dir_folds ,'dir')
            mkdir(dir_folds);
        end
        
%         try
            Arch_file_name = ['Archives/' algo_name{k} '/' func_name{func_num} '_' num2str(func_config.Xdim) 'd_' num2str(algo_config.PopuSize) '.mat'];
            file_name = [ dir_folds func_name{func_num} '_' num2str(func_config.Xdim) 'd_' num2str(algo_config.PopuSize) '.mat'];
            if exist(Arch_file_name,'file') && ~strcmp(result_folder ,  'Archives/')
                copyfile(Arch_file_name , file_name);
                disp(['Find ' file_name ' in Archives']);
            elseif ~exist(file_name,'file')
                alrecord(algo_name{k}, func_name{func_num} , algo_config,func_config, totle_RUN , file_name);
            else
                disp(['  already exist' file_name ' need not to calculate']);
            end
%         catch ErrorInfo
%                 disp(ErrorInfo);  
%         %     disp(ErrorInfo.identifier);  
%             disp(ErrorInfo.message);  
%         %     disp(ErrorInfo.stack);  
%         %     disp(ErrorInfo.cause);  
%             continue;
%         end
    end
end

for func_num = 1:nFunc_num
    % **************************==- analysis -==***************************
    % ********==- draw online figure -==********
    draw_data = cell(nAlgorithm);
    for  k = 1:nAlgorithm
        draw_data{k} = [];
        file_name = [result_folder  algo_name{k} '/' func_name{func_num} '_' num2str(func_config.Xdim) 'd_' num2str(algo_config.PopuSize) '.mat'];
        [ result_series  , run_series , Result_val ]  = analysis( file_name , Run ,totle_RUN);
        Fx(:,k) = result_series;
        draw_data{k} = run_series;
        finalRes( (func_num-1)*nAlgorithm + k , : ) = [Result_val 0];
        % draw
        leng_str{k} = ['\fontsize{15}' algo_name{k}];
    end
    close all;
    draw_lines(draw_data , nAlgorithm, Run , func_num ,leng_str,1,result_folder);

    % ********==- Wilcoxon -==********
    %  Wilcoxon rank-sum test
    for j=1:nAlgorithm     % ttest
        [P,H,STATS] = ranksum(Fx(:,j), Fx(:,nAlgorithm));
        finalRes((func_num-1)*nAlgorithm+j,5) = sign(finalRes(func_num*nAlgorithm,2) - finalRes((func_num-1)*nAlgorithm+j,2))*H;
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
result1 = Sr;
result2 = finalRes;
result3 = Transpose;

save( [result_folder 'result.mat'] , 'finalRes','result1','result2','result3');


display(['Total time ',num2str(cputime - totaltime) ' s'] );

