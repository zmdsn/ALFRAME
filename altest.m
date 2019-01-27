function [  ] = altest( algo , func , is_figure)
% 算法测试函数
% algo 算法名称
% func 测试函数
    addroad;
%     clc; 
    
    %% 选择算法
    algo_name = algo;
    % 设置算法参数
    algo_config.PopuSize = 40;
    algo_config.maxFES = 60000;

    %% 选择测试函数
    func_name = func;
    % func_name = 'mobjf1';
    
    % 设置函数参数
    func_config.Xdim = 30; 

    % 测试一个算法
    [best_fvalue,best_solution,run_series,run_info] = feval(algo_name,func_name,algo_config,func_config);
%     [best_solution,best_fvalue,run_series,run_info] = feval(algo_name,func_name,algo_config);
    if exist('is_figure','var')
        draw_data{1} = run_series;
        draw_lines(draw_data , 1, 1 , 1 ,[algo_name],1);
    end
    disp( ['  ************** ' algo_name ' ***************'] );
    disp(['      Best Value  =  ' num2str(best_fvalue,'%10.5e') ]);
    disp(['      Time Used   =  ' num2str(run_info) 's']);
end

