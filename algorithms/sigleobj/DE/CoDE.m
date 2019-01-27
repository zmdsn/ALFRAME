function [best_fvalue,best_solution,run_series,run_info,user_set] = CoDE(func_name , algo_config,func_config,user_config)

%**************************************************************************************************
%Author: Yong Wang
%Last Edited: July 1, 2010
%Email: ywang@csu.edu.cn
%Reference: Differential Evolution with Composite Trial Vector Generation Strategies
%                                                    and Control Parameters
%                           IEEE Transactions on Evolutionary Computation, Accepted
%
% Adapte Author : Algo
% Email : zmdsn@126.com
% Date : 9/28/2016
%**************************************************************************

    %% ****************==- Initialization settings -==***********************
    format long;
    format compact;
    rng('shuffle'); 
    % set default
    default_set = struct('PopuSize' , 50,...
                         'maxFES'   , 300000  );
                     
    if ~exist('func_name','var') func_name = @sin; end
    if ~exist('algo_config','var') algo_config = []; end
    if ~exist('func_config','var') func_config = []; end
    % get parameters
    [algo_config, func_config, objfunc, popsize, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
    %% ****************==- main body  -==***********************

    lu = [Xmin * ones(1, DIM); Xmax * ones(1, DIM)];
    time = 1;
    outcome = [];
    % The total number of runs
    totalTime = 1;
    %% *****************==- Main body which was provided by the authors -==**********************
    while time <= totalTime

        % Initialize the main population
        p = repmat(lu(1, :), popsize, 1) + rand(popsize, DIM) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));

        % Evaluate the objective function values
        fit = f(p',func_config);
        if(size(fit,1) == 1)
            fit = fit';               
        end
        % Record the number of function evaluations (FES)
        FES = popsize;
        while FES < maxFES
            pTemp = p;
            fitTemp = fit;           
            % uSet: the set of trial vectors
            uSet = zeros(3 * popsize, DIM);
            for i = 1 : popsize
                % the three control parameter settings
                F  = [1.0 1.0 0.8];
                CR = [0.1 0.9 0.2];
                % Uniformly and randomly select one of the control
                % parameter settings for each trial vector generation strategy
                paraIndex = floor(rand(1, 3) * length(F)) + 1;

                % Generate the trail vectors
                u = generator(p, lu, i, F, CR, popsize, DIM, paraIndex);
                uSet(3*i -2 : 3*i, :) = u;
                FES = FES + 3;
            end
            % Evaluate the trial vectors
            fitSet = f(uSet',func_config);
            if(size(fitSet,1) == 1)
                fitSet = fitSet';               
            end
            for i = 1 : popsize
                % Choose the best trial vector from the three trial vectors

                [minVal, minID] = min(fitSet(3*i-2 : 3*i, :));
                bestInd = uSet(3*(i-1) + minID, :);
                bestIndFit = fitSet(3*(i-1) + minID, :);
                % Choose the better one between the trial vector and the
                % target vector
                if fit(i) >= bestIndFit
                    pTemp(i, :) = bestInd;
                    fitTemp(i, :) = bestIndFit;
                end
            end
            p = pTemp;
            fit = fitTemp;
            if FES > maxFES
                FES = maxFES;
            end
            outcome = [outcome; [FES min(fit)]];
        end
        time = time + 1;
    end
    %% ****************==- 算法结果数据处理 -==*********************
    [best_fvalue,minidx]= min(fit);
    best_solution = uSet(minidx, :)';
    run_series = outcome;
    Altime = cputime - mytime;                  % 计算时间
    run_info = [Altime];
end

function u = generator(p, lu, i, F, CR, popsize, n, paraIndex)

    %.... "rand/1/bin" strategy ....%

    % Choose the indices for mutation
    indexSet = 1 : popsize;
    indexSet(i) = [];

    % Choose the first Index
    temp = floor(rand * (popsize - 1)) + 1;
    index(1) = indexSet(temp);
    indexSet(temp) = [];

    % Choose the second index
    temp = floor(rand * (popsize - 2)) + 1;
    index(2) = indexSet(temp);
    indexSet(temp) = [];

    % Choose the third index
    temp = floor(rand * (popsize - 3)) + 1;
    index(3) = indexSet(temp);

    % Mutation
    v1 = p(index(1), :) + F(paraIndex(1)) .* (p(index(2), :) - p(index(3), :));

    % Handle the elements of the mutant vector which violate the boundary
    vioLow = find(v1 < lu(1, :));
    if ~isempty(vioLow)
        v1(1, vioLow) = 2 .* lu(1, vioLow) - v1(1, vioLow);
        vioLowUpper = find(v1(1, vioLow) > lu(2, vioLow));
        if ~isempty(vioLowUpper)
            v1(1, vioLow(vioLowUpper)) = lu(2, vioLow(vioLowUpper));
        end
    end

    vioUpper = find(v1 > lu(2, :));
    if ~isempty(vioUpper)
        v1(1, vioUpper) = 2 .* lu(2, vioUpper) - v1(1, vioUpper);
        vioUpperLow = find(v1(1, vioUpper) < lu(1, vioUpper));
        if ~isempty(vioUpperLow)
            v1(1, vioUpper(vioUpperLow)) = lu(1, vioUpper(vioUpperLow));
        end
    end

    % Binomial crossover
    j_rand = floor(rand * n) + 1;
    t = rand(1, n) < CR(paraIndex(1));
    t(1, j_rand) = 1;
    t_ = 1 - t;
    u(1, :) = t .* v1 + t_ .* p(i, :);

    %... "current to rand/1" strategy ...%

    % The mechanism to choose the indices for mutation is slightly different from that of the classic
    % "current to rand/1", we found that using the following mechanism to choose the indices for
    % mutation can improve the performance to certain degree
    index(1: 3) = floor(rand(1, 3) * popsize) + 1;

    % Mutation
    v2 = p(i, :) + rand * (p(index(1), :) - p(i, :)) + F(paraIndex(2)) .* (p(index(2), :) - p(index(3), :));

    % Handle the elements of the mutant vector which violate the boundary
    vioLow = find(v2 < lu(1, :));
    if ~isempty(vioLow)
        v2(1, vioLow) = 2 .* lu(1, vioLow) - v2(1, vioLow);
        vioLowUpper = find(v2(1, vioLow) > lu(2, vioLow));
        if ~isempty(vioLowUpper)
            v2(1, vioLow(vioLowUpper)) = lu(2, vioLow(vioLowUpper));
        end
    end

    vioUpper = find(v2 > lu(2, :));
    if ~isempty(vioUpper)
        v2(1, vioUpper) = 2 .* lu(2, vioUpper) - v2(1, vioUpper);
        vioUpperLow = find(v2(1, vioUpper) < lu(1, vioUpper));
        if ~isempty(vioUpperLow)
            v2(1, vioUpper(vioUpperLow)) = lu(1, vioUpper(vioUpperLow));
        end
    end

    % Binomial crossover is not used to generate the trial vector under this
    % condition
    u(2, :) = v2;

    %... "rand/2/bin" strategy ...%

    % Choose the indices for mutation
    indexSet = 1 : popsize;
    indexSet(i) = [];

    % Choose the first index
    temp = floor(rand * (popsize - 1)) + 1;
    index(1) = indexSet(temp);
    indexSet(temp) = [];

    % Choose the second index
    temp = floor(rand * (popsize - 2)) + 1;
    index(2) = indexSet(temp);
    indexSet(temp) = [];

    % Choose the third index
    temp = floor(rand * (popsize - 3)) + 1;
    index(3) = indexSet(temp);
    indexSet(temp) = [];

    % Choose the fourth index
    temp = floor(rand * (popsize - 4)) + 1;
    index(4) = indexSet(temp);
    indexSet(temp) = [];

    % Choose the fifth index
    temp = floor(rand * (popsize - 5)) + 1;
    index(5) = indexSet(temp);

    % Mutation
    % The first scaling factor is randomly chosen from 0 to 1
    v3 = p(index(1), :) + rand .* (p(index(2), :) - p(index(3), :)) + F(paraIndex(3)) .* (p(index(4), :) - p(index(5), :));

    % Handle the elements of the mutant vector which violate the boundary
    vioLow = find(v3 < lu(1, :));
    if ~isempty(vioLow)
        v3(1, vioLow) = 2 .* lu(1, vioLow) - v3(1, vioLow);
        vioLowUpper = find(v3(1, vioLow) > lu(2, vioLow));
        if ~isempty(vioLowUpper)
            v3(1, vioLow(vioLowUpper)) = lu(2, vioLow(vioLowUpper));
        end
    end

    vioUpper = find(v3 > lu(2, :));
    if ~isempty(vioUpper)
        v3(1, vioUpper) = 2 .* lu(2, vioUpper) - v3(1, vioUpper);
        vioUpperLow = find(v3(1, vioUpper) < lu(1, vioUpper));
        if ~isempty(vioUpperLow)
            v3(1, vioUpper(vioUpperLow)) = lu(1, vioUpper(vioUpperLow));
        end
    end

    % Binomial crossover
    j_rand = floor(rand * n) + 1;
    t = rand(1, n) < CR(paraIndex(3));
    t(1, j_rand) = 1;
    t_ = 1 - t;
    u(3, :) = t .* v3 + t_ .* p(i, :);
end

