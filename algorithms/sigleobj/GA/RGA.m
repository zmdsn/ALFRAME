%**************************************************************************
% bRGA 
% input Xmin, Xmax 
%
% Author:     algo
% programmed: Apr 18, 2016
% 
%**************************************************************************
function [best_fvalue,best_solution,run_series,run_info,user_set] = bRGA(func_name , algo_config,func_config,user_config)

    %% ****************==- Initialization settings -==***********************
    format long;
    format compact;
    rng('shuffle'); 
    % set default
    default_set = struct('PopuSize' , 100,...
                         'maxFES'   , 300000 );
    % func_config.Xdim = 2;
    if ~exist('func_name','var') func_name = @sin; end
    if ~exist('algo_config','var') algo_config = []; end
    if ~exist('func_config','var') func_config = []; end
    % get parameters
    [algo_config, func_config, objfunc, PopuSize, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
    Xmin = Xmin(1);
    Xmax = Xmax(1);

    % Initial parameters
    outcome = [];    %FES = 0; 
    pCrossover = 0.7;  
    pMuation = 0.1;  % GSA: A Gravitational Search Algorithm
    %myEps = 10^(-4);           %distance of two vectors > myEps,believe are not the same

    % Initialize population
    swarm = objfunc.initialize( PopuSize );

    % Evaluate the objective function values
    f_swarm = objfunc.obj_value( swarm );
    FES = PopuSize;
    [fitval, fitindex] = min(f_swarm);
    globalBestValue = fitval;   
    globalBestSolution = swarm(:,fitindex);
    outcome = [outcome; [FES globalBestValue]];

    newswarm = swarm;   f_newswarm = f_swarm;   
    newswarm1 = swarm;  f_newswarm1 = f_swarm;  
    % Main body of real-coded GA
    while FES < maxFES
        % Selection
        % roulette gambling
        ifs = 1./( f_swarm ) ;
        sumfitness = sum( ifs );
        sumf = ifs ./sumfitness;
        index=[];
        for i=1:PopuSize   %转sizepop次轮盘
            pick=rand;
            while pick==0
                pick=rand;
            end
            for j=1:PopuSize
                pick=pick-sumf(j);
                if pick< 0
                    index=[index j];
                    break;  %寻找落入的区间，此次转轮盘选中了染色体i，注意：在转sizepop次轮盘的过程中，有可能会重复选择某些染色体
                end
            end
        end
        parent = swarm(:,index);
%         f_newswarm = f_swarm(index);
        
        % Crossover
        for mi = 1:(PopuSize)
            randomIndex = randperm(PopuSize);
            % select two solutions for crossover:
            if rand < pCrossover
               ii1 = randomIndex(1); 
               ii2 = randomIndex(2);
               weight = rand(DIM,1)< 0.5; 
%                weight = rand; 
               newswarm(:,mi) = weight.*parent(:,ii1) + (1-weight).*parent(:,ii2); 
            else
               newswarm(:,mi) = parent(:,mi);
            end
        end

        % Mutation:
        for mi = 1 : PopuSize
            for nj =1: DIM
                if rand < pMuation
                   newswarm(nj,mi) = newswarm(nj,mi) + randn;
                end
            end
        end

        % Viloate LB or UB??   % repair operator   
        indexLB = find(newswarm<Xmin); 
        newswarm(indexLB) = min( Xmax, 2*Xmin-newswarm(indexLB) );
        indexUB = find(newswarm1>Xmax);
        newswarm(indexUB) = max( Xmin, 2*Xmax-newswarm(indexUB) );

        % Evaluatation
        f_newswarm = objfunc.obj_value( newswarm );
        FES = FES + PopuSize;

        % Elitism
        [fxopt, g_over] = min(f_newswarm);
        if fxopt <= globalBestValue %even better solution found, update it!!
            globalBestValue = fxopt ;
            globalBestSolution = newswarm(g_over);   
        else % otherwise, replace a random solution if it is not in the swarm
%             [f_worst, g_worst] = max(f_newswarm); % only fitness
            g_worst = randperm(PopuSize,1);
            newswarm(:,g_worst) =  globalBestSolution;
            f_newswarm(g_worst(1)) =  globalBestValue(1); 
        end
        
        % ready for next iteration
        swarm = newswarm;      f_swarm = f_newswarm;
        
       %% 锟斤拷图
%         objfunc.draw_running(swarm);  
        % Update the return best-found-so-far function value
        outcome = [outcome; [FES globalBestValue]];
    end
    % END: function outcome =HypeGA(Xmin, Xmax)
    %% ****************==- 锟姐法锟斤拷锟斤拷锟捷达拷锟斤拷 -==*********************
    best_fvalue = globalBestValue;
    best_solution = globalBestSolution;
    run_series = outcome;
    Altime = cputime - mytime;                  % 锟斤拷锟斤拷时锟斤拷
    run_info = [Altime];
end