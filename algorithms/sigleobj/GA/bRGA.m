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
    default_set = struct('PopuSize' , 5,...
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
    pCrossover = 0.8;  pMuation = 0.1;  % GSA: A Gravitational Search Algorithm
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
        % Crossover
        randomIndex = randperm(PopuSize);
        for mi = 1:2:(PopuSize-1)
            % select two solutions for crossover:
            if rand < pCrossover
               ii1 = randomIndex(mi); 
               ii2 = randomIndex(mi+1);
               weight = rand; 
               %weight = f_swarm(ii2)/(f_swarm(ii1)+f_swarm(ii2)); % maybe better, adaptive
               newswarm(:,ii1) = weight*swarm(:,ii1) + (1-weight)*swarm(:,ii2); 
               newswarm(:,ii2) = (1-weight)*swarm(:,ii1) + weight*swarm(:,ii2); 
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

        % Selection
          %1: Tournament selection
    %      for mi = 1 : PopuSize
    %          tourIndex = randperm(PopuSize);
    %          tourIndex = tourIndex(1:tournamentSize);
    %          [f_newswarm1(mi),in14] =  min(f_newswarm(tourIndex));
    %          newswarm1(mi,:) = newswarm(tourIndex(in14),:);
    %          f_newswarm1(mi) =  f_newswarm(tourIndex(in14));
    %      end
           %2: (lamda + mu)
         f_union = [f_swarm  f_newswarm ];
         [tempValue, tempIndex] = sort(f_union);
         for mi = 1 : PopuSize
            if tempIndex(mi)<=PopuSize 
                newswarm1(:,mi) = swarm(:,tempIndex(mi));
                f_newswarm1(mi) = f_swarm(tempIndex(mi));
            else
                newswarm1(:,mi) = newswarm(: , tempIndex(mi)-PopuSize);
                f_newswarm1(mi) = f_newswarm(tempIndex(mi)-PopuSize);
            end
         end
         
        % Elitism
        [fxopt, g_over] = min(f_newswarm1);
        if fxopt <= globalBestValue %even better solution found, update it!!
            globalBestValue = fxopt ;
            globalBestSolution = newswarm1(:,g_over);   
        else % otherwise, replace a random solution if it is not in the swarm
            Distance = ones(1, PopuSize);
            for mi=1:PopuSize
                Distance(mi) = norm(newswarm1(:,mi)-globalBestSolution);
            end
            [f_worst, g_worst] =  max(f_newswarm1(1:PopuSize)./(Distance+eps));       
            %[f_worst, g_worst] = max(f_swarm); % only fitness
            newswarm1(:,g_worst) =  globalBestSolution;
            f_newswarm1(g_worst(1)) =  globalBestValue(1); 
        end

        
% 
%         B = pinv(swarm)*newswarm1
%         det(B)
        % ready for next iteration
        swarm = newswarm1;      f_swarm = f_newswarm1;   
       %% ��ͼ
%         objfunc.draw_running(swarm);  
        % Update the return best-found-so-far function value
        outcome = [outcome; [FES globalBestValue]];
    end
    % END: function outcome =HypeGA(Xmin, Xmax)
    %% ****************==- �㷨�����ݴ��� -==*********************
    best_fvalue = globalBestValue;
    best_solution = globalBestSolution;
    run_series = outcome;
    Altime = cputime - mytime;                  % ����ʱ��
    run_info = [Altime];
end