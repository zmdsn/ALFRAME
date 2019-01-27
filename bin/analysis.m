function [ result_series , online_series , Result_val ] = analysis( filename , Run ,totle_RUN)
% return the result from saved data
% algo : the name of algorithm
% func : the name of function
    
    run_data = [];
    seled = Run;

    load( filename );
    % select seled result from run
    selidx = randperm(totle_RUN);
    selidx = selidx(1:seled) + 1;
    
    % mean online series    
    online_series(:,1) = run_data(:,1);

    online_series(:,2) = mean(run_data(:,selidx)')';
    
    % the result series 
    result_series = run_data(end,selidx);

    Result_val = [min(result_series)',mean(result_series)',median(result_series)',std(result_series)'];
    
end

