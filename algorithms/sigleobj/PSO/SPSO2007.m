function [best_fvalue,best_solution,run_series,run_info,user_set] = SPSO2007(func_name , algo_config,func_config,user_config)

%**************************************************************************
% Abstract :  
% Standard Particle Swarm Optimisation
%
% Standard PSO 2007
% Developed by: Dr. Mahamed G.H. Omran (omran.m@gust.edu.kw) 
% 
% Adapted Author : Algo
% Email : zmdsn@126.com
% Date : 6/13/2017
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
    [algo_config, func_config, objfunc, PopuSize, maxFES, DIM, Xmin, Xmax, func, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
    D = DIM;          % 获取目标函数类的一些参数
    LB =  objfunc.Xmin(1);
    UB = objfunc.Xmax(1);
    N = PopuSize ;
    FE_max = maxFES;
    %% ****************==- main body  -==***********************

    outcome = [];
    err = 1e-365;
    % Standard PSO 2007
    % Developed by: Dr. Mahamed G.H. Omran (omran.m@gust.edu.kw) 7-May-2011

    % Inputs: -----
    % x : swarm
    % f : swarm fitness
    % FE_max : maximum number of function evaluations 
    % fun : function specifier
    % err : admissible error
    % LB : Lower bound
    % UB : Upper bound
    % opt_f : global optimum for fun

    % Outputs: -----
    % best_x : best solution
    % best_f : least error (best fitness - opt_f)
    % max_FEs : number of function evaluations to reach the admissible error

    % In order to run this program you have to create a function called "func"
    % that has the following structure

    % function [fit] = func(x, fun)
    % where x is a D-dimensional vector representing a particle and fun is a
    % function specifier. Within the body of this function you have to do
    % something like this:
    % 
    % D = length(x);
    % 
    % switch(fun)
    %     case {1}, % Sphere function
    %         fit = sum(x.^2);
    %     case {2}, % Rastrigin
    %         fit = 10*D + sum(x.^2 - 10*cos(2.*pi.*x));
    % end
    % end


    % Setting PSO parameters
    w = 1/(2*log(2));
    c1 = 0.5 + log(2);
    c2 = c1;

    % N : swarm size (should be set to floor(10 + 2*sqrt(D)))
    % D : problem dimension
    x = init( N )';
    % personal best experience
    p_x = zeros(size(x));
    % personal best fitness
    p_f =  func(x',func_config);
    % velocity
    v = zeros(size(x));

    % initialize the personal experience
    % p_x = x;
    % p_f = f;


    % index of the global best particle
    [best_f, g] = min(p_f);

    % initialize velocities
    for i=1:1:N
        r_p = LB + rand(1,D) .* (UB - LB);
        v(i,:) = 0.5*(r_p - x(i,:));
    end

    %PSO_run = zeros(1,t_max); %a PSO run
    FEs = 0;
    max_FEs = 0;
    count = 1;

    % K neighbors for each particle -- based on Clerc description 
    % http://clerc.maurice.free.fr/pso/random_topology.pdf
    % P. 2 (Method 2)
    K = 3;

    while FEs < FE_max

        % In the C version, random permutation is applied here. This is
        % currently not implemented in this code.

            if count > 0  % No improvement in the best solution. So randomize topology

                L = eye(N,N);

                for s = 1:1:N
                    for k=1:1:K
                        rr = 1 + rand*(N - 1);
                        r = floor(rr + 0.5);
                        L(s,r) = 1;
                    end
                end

            end % if


        for i = 1:1:N

            % Find the best informant

            MIN = Inf;

            for s=1:1:N
                if (L(s,i) == 1)
                    if p_f(s) < MIN
                        MIN = p_f(s);
                        g_best = s;
                    end
                end
            end

            % Velocity update equation    
            v(i,:) = w*v(i,:) + c1*rand(1,D).*(p_x(i,:) - x(i,:)) + c2*rand(1,D).*(p_x(g_best,:) - x(i,:));

            % Position update equation
            x(i,:) = x(i,:) + v(i,:);

           %Check for constraint violations
           for j = 1:1:D

               if x(i,j) > UB
                   x(i,j) = UB;
                   v(i,j) = 0;
               end

               if x(i,j) < LB
                   x(i,j) = LB;
                   v(i,j) = 0;
               end

           end %j

           % calculate fitness of particle i
            f(i) = func(x(i,:)',func_config);

            FEs = FEs + 1;

        end %i

        % Update personal best
        for i=1:1:N
            if f(i) <= p_f(i)
                p_x(i,:) = x(i,:);
                p_f(i) = f(i);
            end %if
        end %i

        % Update global best
        [b_f, g] = min(p_f);
        if b_f < best_f
            best_f = b_f;
            
            count = 0;
        else
            count = count + 1;
        end
        outcome = [outcome; [FEs best_f]];
    end %t

    best_x = p_x(g,:);

    if max_FEs == 0
        max_FEs = FEs;
    end
    
    %% ****************==- collating the results -==*********************
    best_fvalue = best_f;
    best_solution = best_x;
    run_series = outcome;
    Altime = cputime - mytime;                  % 计算时间
    run_info = [Altime];
    
end %PSO

