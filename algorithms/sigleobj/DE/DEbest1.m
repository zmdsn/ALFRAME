function [best_fvalue,best_solution,run_series,run_info,user_set] = DErand1(func_name , algo_config,func_config,user_config)

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
                         'maxFES'   , 3000  );
                     
    if ~exist('func_name','var') func_name = @sin; end
    if ~exist('algo_config','var') algo_config = []; end
    if ~exist('func_config','var') func_config = []; end
    % get parameters
    [algo_config, func_config, objfunc, PopuSize, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
    %% ****************==- main body  -==***********************


    FES = 0;
    outcome = [];    
    globalBestValue = Inf;
    swarm = init( PopuSize );  % ��ʼ����Ⱥ,���һ����Ⱥ��ģΪPopuSize����Ⱥ
    n_Strategies = 1;
    f_swarm = f( swarm ,func_config);
    Vswarm = swarm;  f_Vswarm = f_swarm;   
    Uswarm = swarm;  f_Uswarm = f_swarm; 
    % Initial parameters
    outcome = [];    %FES = 0; 
    % F: the scaling factor of DE: different F
    F = normrnd(0.5, 0.3, 1, PopuSize); %0.5 + rand*0.1;
    % CR: the crossover control parameter of DE
    rr = (rand(1, DIM)>0.5);  %different CR with 0.9 or 0.2
    CR = 0.9.*rr + 0.2.*(1-rr); %0.9 + rand * 0.05;
    [fitval, fitindex] = min(f_swarm);
    globalBestValue = fitval;   
    globalBestSolution = swarm(:,fitindex );
    outcome = [outcome; [FES globalBestValue]];
    FES = FES + PopuSize;
    
    % ��ѭ��
    while FES < maxFES
       %%  1 Mutation
        for mi = 1:PopuSize
            choice = 1 + floor(rand*n_Strategies);
            switch choice
                case 1
                    rIn = RandIndex(PopuSize,2);
                    Vswarm(:,mi) = globalBestSolution+F(mi)*(swarm(:,rIn(1))-swarm(:,rIn(2)));
                otherwise
                   display([' no this choice = ' num2str(choice)]);
            end
        end
        %%
        %%2  Binomial Crossover
        for mi = 1:PopuSize
            jRand = floor(rand*DIM) + 1;
            posCross = (rand(1,DIM)<CR);
            posCross(jRand) = 1;
            posCross_ = 1 - posCross;
            Uswarm(:,mi) = posCross'.*Vswarm(:,mi) + posCross_'.*swarm(:,mi);
        end
        %2 Exponential Crossover

        %%
        %3 Viloate LB or UB??   % repair operator   
        indexLB = find(Uswarm<Xmin); 
        Uswarm(indexLB) = min( Xmax, 2*Xmin-Uswarm(indexLB) );
        indexUB = find(Uswarm>Xmax);
        Uswarm(indexUB) = max( Xmin, 2*Xmax-Uswarm(indexUB) );

        %4 Evaluatation
        f_Uswarm = f( Uswarm ,func_config); % ����Ŀ�꺯��ֵ,ע�������Ǻ���ֵ������Ӧ��ֵ        [fitval, fitindex] = min(f_swarm);   % �����Сֵ,�˴�����Ҫʹ�õ����� sort()
        FES = FES + PopuSize;     
       %%
        %5 Selection: store into "Vswarm" temporarily for Elitism
         for mi = 1 : PopuSize
            if f_Uswarm(mi)<= f_swarm(mi) 
                Vswarm(:,mi) = Uswarm(:,mi);
                f_Vswarm(mi) = f_Uswarm(mi);
            else
                Vswarm(:,mi) = swarm(:,mi);
                f_Vswarm(mi) = f_swarm(mi);
            end
         end    
        %%
        [fitval, fitindex] = min(f_swarm);
        globalBestValue = fitval(1,1);           
        globalBestSolution = swarm(:,fitindex(1)); 
        %6 Elitism
        [fxopt, g_over] = min(f_Vswarm); 
        if fxopt <= globalBestValue %even better solution found, update it!!
            globalBestValue = fxopt;   
            globalBestSolution = Vswarm(:,g_over );   
        else % otherwise, replace a random solution if it is not in the swarm
            Distance = ones(1, PopuSize);
            for mi=1:PopuSize
                Distance(mi) = norm(Vswarm(:,mi)-globalBestSolution);
            end
            [f_worst, g_worst] =  max(f_Vswarm(1:PopuSize)./(Distance+eps));       
            Vswarm(:,g_worst) =  globalBestSolution;
            f_Vswarm(g_worst) =  globalBestValue; 
        end
        swarm = Vswarm;      f_swarm = f_Vswarm;   
        % ��ͼ
%         hold off;
%         scatter(swarm(1,:),swarm(2,:));
%         hold on;
%         objfunc.fun_figue(3);   
%         pause(.1);
%         
        outcome = [outcome; [FES globalBestValue]]; % �㷨���Ϊ���еľ���
    end
    %% ****************==- collating the results -==*********************
    best_fvalue = globalBestValue;
    best_solution = globalBestSolution;
    run_series = outcome;
    Altime = cputime - mytime;                  % ����ʱ��
    run_info = [Altime];
end

function randindex = RandIndex(pSize,k)
%size: input popusize;  k: k different subindex are required
%randindex: k different subindex are returned in this array
randpert = randperm(pSize);
randindex = randpert(1:k);
end