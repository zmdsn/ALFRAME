function [best_fvalue,best_solution,run_series,run_info,user_set] = TPSO(func_name , algo_config,func_config,user_config)

%**************************************************************************
% Abstract :  
    % Kennedy J, Eberhart R. Particle swarm optimization, in Neural
    % Networks, 1995.In: Proceedings, IEEE international conference on; 1995. p. 1942�C1948.
    % use 
%
% Author : Algo
% Email : zmdsn@126.com
% Date : 9/28/2016
%**************************************************************************

    %% ****************==- Initialization settings -==***********************
    format long;
    format compact;
    rng('shuffle'); 
    % set default
    default_set = struct('PopuSize' , 10,...
                         'maxFES'   , 300000);
                     
    if ~exist('func_name','var') func_name = @sin; end
    if ~exist('algo_config','var') algo_config = []; end
    if ~exist('func_config','var') func_config = []; end
    func_config.Xdim = 10;
    % get parameters
    [algo_config, func_config, objfunc, PopuSize, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
 

    %% ****************==- main body  -==***********************
    % **********==- Initialization algorithm parameter -==*******
    FES = 0;
    outcome = [];

    swarm = init( PopuSize );  % ��ʼ����Ⱥ,���һ����Ⱥ��ģΪPopuSize����Ⱥ
    bestL = Xmin + (Xmax - Xmin).*rand(DIM,PopuSize);
    v = zeros(DIM,PopuSize);   %generate initial searching directions
    newPAR = swarm;     %newPAR:  temporary population
    w = 1/(2*log(2));
    c1 = 0.5 + log(2);
    c2 = c1;
    globalBestValue = Inf;
    f_swarm = zeros(1,PopuSize);
    SS = {};
    DD = {};
    ee = 1;
    %% *********************==- main loop -==************************
    while FES < maxFES
        if FES > PopuSize
            old_f_swarm = f_swarm;
        end
        
        f_swarm = f( swarm , func_config); % ����Ŀ�꺯��ֵ,ע�������Ǻ���ֵ������Ӧ��ֵ
        [fitval, fitindex] = min(f_swarm);
        if FES > PopuSize
            iii = find(f_swarm <= old_f_swarm);  % updating local best solution
            bestL(:,iii) = swarm(:,iii);
        end
        if (fitval(1,1) < globalBestValue)
            globalBestValue = fitval(1,1);            
            globalBestSolution = swarm(:,fitindex(1));  
        end
        
        FES = FES + PopuSize;                       
        outcome = [outcome; [FES globalBestValue]]; 

        if FES >= maxFES
            best_solution = globalBestSolution;
            best_fvalue = globalBestValue;
            break; 
        end
            %% draw
%         objfunc.draw_running(swarm,'PSO');
        
        % updating the velocity and position
        for i=1:PopuSize
              v(:,i) = w*v(:,i) + c1*rand(DIM,1).*(bestL(:,i)-swarm(:,i)) + c2*rand(DIM,1).*(globalBestSolution(:,1)-swarm(:,i));
              newPAR(:,i) = swarm(:,i) + v(:,i);
        end
        
        % Out of bounds      
        indexLB = find( newPAR < Xmin );
        newPAR(indexLB) = min( Xmax(indexLB), 2*Xmin(indexLB) - newPAR(indexLB) );
        indexUB = find( newPAR > Xmax );
        newPAR(indexUB) = max( Xmin(indexUB), 2*Xmax(indexUB) - newPAR(indexUB) );
        
        DD{ee} = newPAR;
        B = pinv(swarm)*newPAR;
        det(B);
        SS{ee} = B;
        ee = ee+1;
        swarm = newPAR;
        
%         objfunc.draw_running( swarm );

    end
    save('SS.mat','SS','DD')
    %% ****************==- collating the results -==*********************
    run_series = outcome;
    Altime = cputime - mytime;                  % ����ʱ��
    run_info = [Altime];
end