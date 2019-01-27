function [best_fvalue,best_solution,run_series,run_info,user_set] = DCGA(func_name , algo_config,func_config,user_config)
%**************************************************************************
% bRGA 
% input Xmin, Xmax 
%
% Author:     algo
% programmed: Apr 18, 2016
% 
%**************************************************************************

    format long;
    format compact;
    rng('shuffle'); 
    % set default
    default_set = struct('PopuSize' , 100,...
                         'maxFES'   , 300000  );
    if ~exist('user_config','var')  user_config.M = [0.56 0.10 0.34]; end
    if ~exist('func_name','var') func_name = @sin; end
    if ~exist('algo_config','var') algo_config = []; end
    if ~exist('func_config','var') func_config = []; end
    % get parameters
    [algo_config, func_config, objfunc, PopuSize, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
    M = user_config.M;
    
    % Initial parameters
    outcome = [];    %FES = 0; 
     % GSA: A Gravitational Search Algorithm
    %myEps = 10^(-4);           %distance of two vectors > myEps,believe are not the same

    % Initialize population
    swarm = init( PopuSize );
%     min(abs(swarm'))
    
% bar(swarm)
    % Evaluate the objective function values
    f_swarm = f( swarm , func_config);
    FES = PopuSize;
    [fitval, fitindex] = min(f_swarm);
    globalBestValue = fitval;   
    globalBestSolution = swarm(:,fitindex);
    outcome = [outcome; [FES globalBestValue]];
    newswarm = swarm;   f_newswarm = f_swarm;   
    warning('off')
    
    
    bestL = swarm;
    f_bestL = f_swarm;
    globalBestSolution = swarm(:,1);
    globalBestValue = f_swarm(1);
    G = swarm;
    KK = 1;
    % PSO
    w = 1/(2*log(2));
    c1 = 0.5 + log(2);
    c2 = c1;
    v = zeros(DIM,PopuSize);
    
    
    
    % Main body of real-coded GA
    while FES < maxFES
        r = rand;
        % generate W
        if r < M(1)
            % generate by GA
            W = getWGA(swarm,f_swarm,PopuSize,DIM);
        elseif r < M(2)
            % generate by PSO
            [W,v] = getWPSO(swarm,v,bestL,globalBestSolution,PopuSize,DIM);
        else
            % generate by DE
            W = getWDE(swarm,PopuSize,DIM,globalBestSolution)    ;
        end
%             % generate by self
%             KK = 0.01;
%             G = repmat(globalBestSolution,1,PopuSize).*((rand(DIM,PopuSize)-0.5)*KK + 1);
%             WKK = swarm\G;
%             
% %         W = M(1)*WGA + M(2)*WPSO + M(3)*WDE + M(4)*WKK ;
%         W = M(1)*WGA + M(2)*WPSO + M(3)*WDE ;
        newswarm = swarm*W;

        
        % evaluate newswarm
        f_newswarm = f( newswarm , func_config); % 

        % elitist choice and updating the local best and global best
        idx = find( f_newswarm <= f_swarm );
        swarm(:,idx) = newswarm(:,idx);
        f_swarm(idx) = f_newswarm(idx);
        
%         objfunc.draw_running( swarm );

        [fitval, fitindex] = min(f_swarm);        
        % updating local best solution
        iii = find(f_swarm <= f_bestL);  
        bestL(:,iii) = swarm(:,iii);
        f_bestL(iii) = f_swarm(iii);
        
        % updating global best solution
        if (fitval(1,1) < globalBestValue)
            globalBestValue = fitval(1,1);            
            globalBestSolution = swarm(:,fitindex(1));  
        end
        
        if FES >= maxFES
            best_solution = globalBestSolution;
            best_fvalue = globalBestValue;
            break; 
        end
        
        FES = FES + PopuSize;
        [globalBestValue,idy] = min(f_swarm);
        globalBestSolution = swarm(:,idy);
        outcome = [outcome; [FES globalBestValue]];

    end
%     a
    % END: function outcome =HypeGA(Xmin, Xmax)

    %% ****************==- collating the results -==*********************
    best_fvalue = globalBestValue;
    best_solution = globalBestSolution;
    run_series = outcome;
    Altime = cputime - mytime;                  
    run_info = [Altime];
end










