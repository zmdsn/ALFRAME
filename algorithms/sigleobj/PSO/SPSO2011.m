function [best_fvalue,best_solution,run_series,run_info,user_set] = SPSO2011(func_name , algo_config,func_config,user_config)

%**************************************************************************
% Abstract :  
% Standard Particle Swarm Optimisation
% 
% Standard PSO 2011
% Developed by: Dr. Mahamed G.H. Omran (omran.m@gust.edu.kw) 7-May-2011
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
    LB =  objfunc.Xmin;
    UB = objfunc.Xmax;
    N = PopuSize ;
    FE_max = maxFES;
    %% ****************==- main body  -==***********************

    
    outcome = [];
    % D = objfunc.Xdim;          % 获取目标函数类的一些参数
    % LB =  objfunc.Xmin;
    % UB = objfunc.Xmax;
    err = 1e-365;
    % x swarm
    % f swarm fitness
    % N swarm size
    % fun -- function specifier
    % D -- Problem dimension
    % LB -- Lower bound
    % UB -- Upper bound
    % opt_f -- global optimum for fun

    % Compatibility with the original C version
    % 
    % 1) RNG and numerical instability.
    % 
    % We (Clerc and Omran) conducted several experiments on some benchmark functions. 
    % For several functions the results are very similar to the ones of 
    % the C version while for one function (CEC F6) the SR is really different:
    % 
    % F6 results (C code):
    % a) No normalization
    % Avg. fitness = 5.69e+01 (2.06e+02)
    % SR= 49.20 %
    % 
    % b) Normalization
    % Avg. fitness = 6.44e+01 (1.57e+02)
    % SR= 37.2 %
    % 
    % F6 results (Matlab code)
    % a) No normalization
    % 
    % Avg. fitness = 6.05e+01(1.58e+02) 
    % SR = 0%
    % 
    % b) Normalization
    % 
    % Avg. fitness = 5.12e+01(1.37e+02) 
    % SR = 0%
    % 
    % We suspect that there is a problem with the Matlab RNG and/or numerical 
    % instability (we implemented a simple RNG using C and Matlab and run our 
    % programs and still we got different results).
    % 
    % 2) Normalization
    % 
    % It is recommended that you use this option (i.e. randomize = 1) when 
    % the search space in not a hypercube. 
    % If the search space is a hypercube, It is better not normalize 
    % (there is a small difference between the position without any normalisation and the de-normalised one.).
    % 
    % 3) Random permutation
    % 
    % The random permutation of the numbering of the particles before 
    % each step is not included in the Matlab version(usually, it does not 
    % make a great difference in the C version).


        w = 1/(2*log(2));
        c1 = 0.5 + log(2);
        c2 = c1;

        x=zeros(N,D);

        % velocity
        v = zeros(size(x));

        for i=1:1:N
            for j=1:1:D
                xMin=LB(j);
                xMax=UB(j); 
                x(i,j) = alea(xMin,xMax);
                v(i,j) = alea(xMin-x(i,j), xMax-x(i,j));  
            end
        end

        f = zeros(1,N);

        % calculate fitness
        for i = 1:1:N
            f(i) = func(x(i,:)',func_config);
        end


        % initialize the personal experience
        p_x = x;
        p_f = f;


        % index of the global best particle
        [best_f, g] = min(p_f);

        %PSO_run = zeros(1,t_max); %a PSO run
        FEs=N; % N initialisations.  *****FEs = 0;
        %max_FEs = 0;

        count = 1;

        % K neighbors for each particle -- based on Clerc description 
        % http://clerc.maurice.free.fr/pso/random_topology.pdf
        % P. 2 (Method 2)
        K = 3;

        p=1-power(1-1/N,K); % Probability to be an informant

        stop=0;

        while stop<1
            % In the C version, random permutation is applied here. This is
            % currently not implemented in this code.

                if count > 0  % No improvement in the best solution. So randomize topology

              %      L = eye(N,N); % Matlab function, but does not exist in FreeMat
                    L=zeros(N,N);
                    for s=1:1:N
                        L(s,s)=1;
                    end

                    for s = 1:1:N % Each particle (column) informs at most K other at random   
                        for r=1:1:N
                        if (r~=s)
                        if (alea(0,1)<p)
                            L(s,r) = 1;
                        end
                        end
                        end
                    end

                end % if


            for i = 1:1:N  % For each particle (line) ..

                %  ...find the best informant g
                        MIN = Inf;  
                for s=1:1:N
                    if (L(s,i) == 1)
                        if p_f(s) < MIN
                            MIN = p_f(s);
                            g_best = s;
                        end
                    end
                end

                % define a point p' on x-p, beyond p
                p_x_p = x(i,:) + c1*(p_x(i,:) - x(i,:));
                % ... define a point g' on x-g, beyond g
                p_x_l = x(i,:) + c2*(p_x(g_best,:) - x(i,:));

                if (g_best == i) % If the best informant is the particle itself, define the gravity center G as the middle of x-p' 
                    G = 0.5*(x(i,:) + p_x_p);
                else % Usual  way to define G
                    sw = 1/3;
                    G = sw*(x(i,:) + p_x_p + p_x_l);
                end


                rad = norm(G - x(i,:)); % radius = Euclidean norm of x-G       
                x_p = alea_sphere(D,rad)+ G; % Generate a random point in the hyper-sphere around G (uniform distribution)       
                v(i,:) = w*v(i,:) + x_p - x(i,:); % Update the velocity = w*v(t) + (G-x(t)) + random_vector 
                                                                                            % The result is v(t+1)
                x(i,:) = x(i,:) + v(i,:); % Apply the new velocity to the current position. The result is x(t+1)

        %Check for constraint violations
        for j = 1:1:D
            xMin=LB(j);
            xMax=UB(j); 
            if x(i,j) > xMax
                x(i,j) = xMax;
                v(i,j) = -0.5*v(i,j); % variant: 0
            end

            if x(i,j) < xMin
                x(i,j) = xMin;
                v(i,j) = -0.5*v(i,j); % variant: 0
            end   
        end %j

                f(i) = func(x(i,:)',func_config);

                FEs = FEs + 1;
                if (FEs>=FE_max) % Too many FE
                    break; 
                end     
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
                count = count + 1; % If no improvement, the topology will be initialised for the next iteration
            end
            outcome = [outcome; [FEs best_f]];

            if (min(p_f) < err) || (FEs>=FE_max)
              stop=1;  
            end


        end %t


        best_x = p_x(g,:);

        max_FEs = FEs;
    %% ****************==- collating the results -==*********************
    best_fvalue = best_f;
    best_solution = best_x;
    run_series = outcome;
    Altime = cputime - mytime;                  % 计算时间
    run_info = [Altime];
end %SPSO2011


function r= alea (a,b)
% Random real  number between a and b
    
    r= a + rand*(b - a); 
 
end


function [x] = alea_sphere(D, radius)
%*  ******* Random point in a hypersphere ********
% Put  a random point inside the hypersphere S(0,radius) (center 0, radius 1).
% Uniform distribution 

% Developed by: Maurice Clerc (May 2011)
    
x = zeros(1,D); %**

% --------- Step 1. Direction

l = 0; 

for j=1:1:D
   % x(j) = randn; %?
    x(j)=alea_normal(0,1); %**
    l = l + x(j)*x(j);
end
    
l=sqrt(l);

% --------- Step 2. Random radius

r=alea(0,1);

x= r*radius*x/l;
    
end

function y1= alea_normal (mean, std_dev) 
    % Use the polar form of the Box-Muller transformation to obtain a pseudo
    % random number from a Gaussian distribution
    
    % Developed by: Maurice Clerc (May 2011)
    
    w=2;
    
    while w>=1    
        x1 = 2.0 * alea (0, 1) - 1.0;
        x2 = 2.0 * alea (0, 1) - 1.0;
        w = x1 * x1 + x2 * x2;     
    end
    
    w = sqrt (-2.0 * log (w) / w);
    y1 = x1 * w;
    
    if alea(0,1)<0.5
        y1=-y1; 
    end
    
    y1 = y1 * std_dev + mean;
    
end
    