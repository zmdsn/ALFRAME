function [time_used] = alcomp(Run,algo_name,algo_config,func_name,func_config)
totaltime = cputime;   
% ****************==- 框架参数设置 -==*********************
nAlgorithm = length(algo_name);
nFunc_num = length(func_name);
result_folder = ['results/pop_' num2str(algo_config.PopuSize) '_dim_' num2str(func_config.Xdim) '/'];

finalRes = zeros(nFunc_num*nAlgorithm,5);
Result_val = zeros(5);
Algorithm = cell(nAlgorithm,Run);
check_frame; % 框架设置检测
for func_num = 1:nFunc_num
    %算法运行
    fxopt = zeros(1,nAlgorithm); % 存储n个算法一次迭代的结果
    Fx = zeros(Run,nAlgorithm);  % 存储n个算法 nFunc_num 次迭代的结果
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
    
    % ****************==- 绘制函数图像 -==*********************
    % 画图draw_data为Run 次运行的结果总和
    draw_data = cell(nAlgorithm);
    for  k = 1:nAlgorithm
      draw_data{k} = [];
      draw_data{k}(:,1) = Algorithm{k}{1}(:,1);
      draw_data{k}(:,2) = Algorithm{k}{1}(:,2);
      for j = 2 : Run
        draw_data{k}(:,2) = draw_data{k}(:,2) + Algorithm{k}{j}(:,2);
      end
      % 图例
      leng_str{k} = ['\fontsize{15}' algo_name{k}];
    end
    close all;
    draw_lines(draw_data , nAlgorithm, Run , func_num ,leng_str,1,result_folder);
    
    
    % ****************==- 结果数据分析 -==*********************
    Result_val = [min(Fx)',mean(Fx)',median(Fx)',std(Fx)'];
    for j=1:nAlgorithm     % ttest
        [H,P,CI,STATS] = ttest2(Fx(:,j), Fx(:,nAlgorithm),0.05,'both');
        finalRes((func_num-1)*nAlgorithm+j,:) = [Result_val(j,:) STATS.tstat];
    end
    result_xls = [result_folder 'Result.xls'];
    xlswrite(result_xls,finalRes);
end  %endfor func_name

% write to tex.xls
Tex_xls = [result_folder 'Tex.xls'];
tex_list = cell(nFunc_num,1);
for i=1:nFunc_num
    test_func = feval(func_name{i});
    tex_list{i} = test_func.Tex;   
end
xlswrite(Tex_xls,tex_list);

format short;
time_used  = cputime - totaltime;
