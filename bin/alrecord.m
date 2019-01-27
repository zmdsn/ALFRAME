function [  ] = alrecord(algo_name, func_name , algo_config,func_config, run , file_name)
% �㷨���Ժ���
% algo �㷨���?
% func ���Ժ���

%     addroad;

    tic
    if ~exist('run','var') 
        run = 30;
    end
    
    solutions = ones(func_config.Xdim,run);
    [ best_fvalue , best_solution , run_series , run_info] = feval(algo_name,func_name,algo_config,func_config);
    disp( ['  ************** ' func_name ' by ' algo_name  ' Run ' num2str(run) ' ***************'] );
    solutions(:,1) =  best_solution;
    lend = size(run_series,1);
    run_data = zeros(lend,run+1);
    run_data(:,1:2) = run_series;
    t=toc;
    disp(['      Time Used   =  ' num2str(t) 's, Estimated time = ' num2str(t*(run/4+1)) 's' ]);    

    parfor ii  = 2:run
        % main
        [ best_fvalue , best_solution , run_series , run_info] = feval(algo_name,func_name,algo_config,func_config);
%         size(best_solution)
        solutions(:,ii) =  best_solution;
        diff = lend - size(run_series(:,2),1);
%         diff
%         size(run_series(:,2),1)
        if diff > 0
            run_data(:,ii+1) = [run_series(:,2) ; ones(diff,1)*run_series(end,2)];
        else
            run_data(:,ii+1) = run_series(1:lend,2);
        end
    end
    save(file_name,'run_data','solutions');
    
    if exist('is_figure','var')
        draw_data{1} = run_data(:,1:2);
        draw_lines(draw_data , 1, 1 , 1 ,[algo_name],1);
    end
    t=toc;
    disp(['      Time Used   =  ' num2str(t) 's, mean = ' num2str(t/run*4) 's' ]);    
end

