    clear all
    load('CLPSOf1.mat')
    
    default_set = struct('PopuSize' , 100,...
                         'maxFES'   , 300000  );
    if ~exist('user_config','var')  user_config.M = [0.1 0.25 0.65]; end
    if ~exist('func_name','var') func_name = 'f1'; end
    if ~exist('algo_config','var') algo_config = []; end
    if ~exist('func_config','var') func_config = []; end
    % get parameters
    [algo_config, func_config, objfunc, PopuSize, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
     
    swarm = DD{1};
    f_swarm = f( swarm , func_config);
    warning('off')
    
    bestL = swarm;
    f_bestL = f_swarm;
    globalBestSolution = swarm(:,1);
    globalBestValue = f_swarm(1);
    v = zeros(DIM,PopuSize);
  
    
    WGA = getWGA(swarm,f_swarm,PopuSize,DIM);

    % generate by PSO
    [WPSO,v] = getWPSO(swarm,v,bestL,globalBestSolution,PopuSize,DIM);
    
    % generate by DE
    WDE = getWDE(swarm,PopuSize,DIM,globalBestSolution);
    
    W1 = DD{1}\DD{2}; 
    
    
    
    
    
    