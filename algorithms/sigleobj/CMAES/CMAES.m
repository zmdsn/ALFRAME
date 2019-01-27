function [best_fvalue,best_solution,run_series,run_info,user_set] = CMAES(func_name , algo_config,func_config,user_config)

%**************************************************************************************************
% Reference:  Hansen N, Ostermeier A. Completely derandomized self-adaptation in evolution
%                     strategies.[J]. Evolutionary Computation, 2001, 9(2):159.
% author : algo
% Email : zmdsn@126.com
% date : 6/13/2017
%**************************************************************************************************

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
    [algo_config, func_config, objfunc, PopuSize, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
    Xmin = Xmin(1);
    Xmax = Xmax(1);

    
    %% ****************==- main body  -==***********************
        
    % Record the best results
    outcome = [];  
    % Main body which was provided by the authors
    time = 1;
    % The total number of runs
    totalTime = 1;
    FES = 0;
    lu = [Xmin * ones(1, DIM); Xmax * ones(1, DIM)];
%     while time <= totalTime
    while FES < maxFES
        % Initialize the main population

        N = DIM;
        maxeval = maxFES; % stop criteria
         % object parameter start point
        xmeanw = (lu(1, :) + rand(1, N) .* (lu(2, :) - lu(1, :)))';
        % initial step size, minimal step size
        sigma = 0.25;  minsigma = 1e-15;  maxsigma = max(lu(2, :)-lu(1, :))/sqrt(DIM);
        flginiphase = 1;  % Initial phase

        % Parameter setting: selection,
        lambda = 4 + floor(3 * log(N)); mu = floor(lambda/2);

        arweights = log((lambda + 1)/2) - log(1:mu)'; % muXone array for weighted recomb.
        % lambda = 10; mu = 2; arweights = ones(mu, 1); % uncomment for (2_I, 10)-ES
        % parameter setting: adaptation
        cc = 4/(N + 4); ccov = 2/(N + 2^0.5)^2;
        cs = 4/(N + 4); damp = (1 - min(0.7, N * lambda/maxeval)) / cs + 1;

        % Initialize dynamic strategy parameters and constants
        B = eye(N); D = eye(N); BD = B * D; C = BD * transpose(BD);
        pc = zeros(N, 1); ps = zeros(N, 1);
        cw = sum(arweights) / norm(arweights); chiN = N^0.5 * (1 - 1/(4 * N) + 1/(21 * N^2));

        % Generation loop
        %     disp(['  (' num2str(mu) ', ' num2str(lambda) ')-CMA-ES (w = [' num2str(arweights', '%5.2f') '])' ]);
        counteval = 0; flgstop = 0;

        % Boundary
        lb = (ones(lambda, 1) * lu(1, :))';
        ub = (ones(lambda, 1) * lu(2, :))';


        while 1 % if flgstop break;
            
            
            % Set stop flag
            if counteval >= maxeval
                flgstop = 1;
            end

            if flgstop
                break;
            end

            % Generate and evaluate lambda offspring
            arz = randn(N, lambda);
            arx = xmeanw * ones(1, lambda) + sigma * (BD * arz);
            x_ = xmeanw * ones(1, lambda);

            % You may handle constraints here and now.
            % You may either resample columns of arz and/or multiply
            % them with a factor < 1 (the latter will result in a
            % decreased overall step size) and recalculate arx
            % accordingly. Do not change arx or arz in any other
            % way. You may also use an altered arx only for the
            % evaluation of the fitness function, if you leave
            % arx and arz unchanged for the algorithm.

            % Handle the elements of the variable which violate the boundary
            I = find(arx > ub);
            arx(I) = 2 * ub(I) - arx(I);
            aa = find(arx(I) < lb(I));
            arx(I(aa)) = lb(I(aa));
            I = find(arx < lb);
            arx(I) = 2 * lb(I) - arx(I);
            aa = find(arx(I) > ub(I));
            arx(I(aa)) = ub(I(aa));

            U = arx';
            
            arfitness =  f( U' , func_config );
            FES = FES + size(U,1);
            counteval = counteval + lambda;
            % Sort by fitness and compute weighted mean in xmeanw
            [arfitness, arindex] = sort(arfitness); % minimization
            xold = xmeanw; % for speed up of Eq. (14)
            xmeanw = arx(:, arindex(1:mu)) * arweights/sum(arweights);
            zmeanw = arz(:, arindex(1:mu)) * arweights/sum(arweights);

            % Adapt covariance matrix
            pc = (1-cc) * pc + (sqrt(cc * (2-cc)) * cw/sigma) * (xmeanw-xold); % Eq. (14)
            if ~flginiphase % do not adapt in the initial phase
                C = (1-ccov) * C + ccov * pc * transpose(pc);           % Eq. (15)
            end
            % adapt sigma
            ps = (1-cs) * ps + (sqrt(cs * (2-cs)) * cw) * (B * zmeanw);      % Eq. (16)
            sigma = sigma * exp((norm(ps)-chiN)/chiN/damp);        % Eq. (17)

            % Update B and D from C
            if mod(counteval/lambda, 1/ccov/N/5) < 1
                C = triu(C) + transpose(triu(C, 1)); % enforce symmetry
                [B, D] = eig(C);
                % limit condition of C to 1e14 + 1
                if max(diag(D)) > 1e14 * min(diag(D))
                    tmp = max(diag(D))/1e14 - min(diag(D));
                    C = C + tmp * eye(N); D = D + tmp * eye(N);
                end
                D = diag(sqrt(diag(D))); % D contains standard deviations now
                BD = B * D; % for speed up only
            end % if mod

            % Adjust minimal step size
            if sigma * min(diag(D)) < minsigma ...
                    | arfitness(1) == arfitness(min(mu + 1, lambda)) ...
                    | xmeanw == xmeanw ...
                    + 0.2 * sigma * BD(:, 1 + floor(mod(counteval/lambda, N)))
                sigma = 1.4 * sigma;

                % flgstop = 1;
            end
            if sigma > maxsigma
                sigma = maxsigma;
            end

            % Test for end of initial phase
            if flginiphase & counteval/lambda > 2/cs
                if (norm(ps)-chiN)/chiN < 0.05  % step size is not much too small
                    flginiphase = 0;
                end
            end
            
            if FES > maxFES
            	break;
            end
            if mod(FES,size(U,1)*10) == 0
                outcome = [outcome; [counteval arfitness(1)]];
            end
        end % while, end generation loop

        %     disp([num2str(counteval) ': ' num2str(arfitness(1))]);
        %     if exist('sfile', 'var')
        %       disp(['Results saved in ' sfile]);
        %     end

        xmin = arx(:, arindex(1)); % return best point of last generation
        %outcome = [outcome arfitness(1)];
        time = time + 1;

    end
    
    %% ****************==- collating the results -==*********************
    best_fvalue = arfitness(1);
    best_solution = U(arindex(1),:);
    run_series = outcome;
    Altime = cputime - mytime;                  
    run_info = [Altime];
    
     %% ****************==- some useless code  -==*********************
    
%     sort(outcome)
%     mean(outcome)
%     std(outcome)
% end   toc;
% Choose the problems to be tested. Please note that for test functions F7
% and F25, the global optima are out of the initialization range. For these
% two test functions, we do not need to judge whether the variable violates
% the boundaries during the evolution after the initialization.
% problemSet = 1 :  25;
% for problemIndex = func_num
% 
%     problem = problemSet(problemIndex)
% 
%     % Define the dimension of the problem
%     n = 30;
% 
%     switch problem
% 
%         case 1
% 
%             % lu: define the upper and lower bounds of the variables
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             % Load the data for this test function
%             load sphere_func_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 2
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load schwefel_102_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 3
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load high_cond_elliptic_rot_data
%             A = []; a = []; alpha = []; b = [];
% 
%             if n == 2, load elliptic_M_D2,
%             elseif n == 10, load elliptic_M_D10,
%             elseif n == 30, load elliptic_M_D30,
%             elseif n == 50, load elliptic_M_D50,
%             end
% 
%         case 4
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load schwefel_102_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 5
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load schwefel_206_data
%             M = []; a = []; alpha = []; b = [];
% 
%         case 6
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load rosenbrock_func_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 7
% 
%             lu = [0 * ones(1, n); 600 * ones(1, n)];
%             load griewank_func_data
%             A = []; a = []; alpha = []; b = [];
% 
%             c = 3;
%             if n == 2, load griewank_M_D2,
%             elseif n == 10, load griewank_M_D10,
%             elseif n == 30, load griewank_M_D30,
%             elseif n == 50, load griewank_M_D50,
%             end
% 
%         case 8
% 
%             lu = [-32 * ones(1, n); 32 * ones(1, n)];
%             load ackley_func_data
%             A = []; a = []; alpha = []; b = [];
% 
%             if n == 2, load ackley_M_D2,
%             elseif n == 10, load ackley_M_D10,
%             elseif n == 30, load ackley_M_D30,
%             elseif n == 50, load ackley_M_D50,
%             end
% 
%         case 9
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load rastrigin_func_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 10
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load rastrigin_func_data
%             A = []; a = []; alpha = []; b = [];
%             if n == 2, load rastrigin_M_D2,
%             elseif n == 10, load rastrigin_M_D10,
%             elseif n == 30, load rastrigin_M_D30,
%             elseif n == 50, load rastrigin_M_D50,
%             end
% 
%         case 11
% 
%             lu = [-0.5 * ones(1, n); 0.5 * ones(1, n)];
%             load weierstrass_data
%             A = []; a = []; alpha = []; b = [];
%             if n == 2, load weierstrass_M_D2, ,
%             elseif n == 10, load weierstrass_M_D10,
%             elseif n == 30, load weierstrass_M_D30,
%             elseif n == 50, load weierstrass_M_D50,
%             end
% 
%         case 12
% 
%             lu = [-pi * ones(1, n); pi * ones(1, n)];
%             load schwefel_213_data
%             A = []; M = []; o = [];
% 
%         case 13
% 
%             lu = [-3 * ones(1, n); 1 * ones(1, n)];
%             load EF8F2_func_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 14
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load E_ScafferF6_func_data
%             if n == 2, load E_ScafferF6_M_D2, ,
%             elseif n == 10, load E_ScafferF6_M_D10,
%             elseif n == 30, load E_ScafferF6_M_D30,
%             elseif n == 50, load E_ScafferF6_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 15
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func1_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 16
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func1_data
%             if n == 2, load hybrid_func1_M_D2,
%             elseif n == 10, load hybrid_func1_M_D10,
%             elseif n == 30, load hybrid_func1_M_D30,
%             elseif n == 50, load hybrid_func1_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 17
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func1_data
%             if n == 2, load hybrid_func1_M_D2,
%             elseif n == 10, load hybrid_func1_M_D10,
%             elseif n == 30, load hybrid_func1_M_D30,
%             elseif n == 50, load hybrid_func1_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 18
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func2_data
%             if n == 2, load hybrid_func2_M_D2,
%             elseif n == 10, load hybrid_func2_M_D10,
%             elseif n == 30, load hybrid_func2_M_D30,
%             elseif n == 50, load hybrid_func2_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 19
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func2_data
%             if n == 2, load hybrid_func2_M_D2,
%             elseif n == 10, load hybrid_func2_M_D10,
%             elseif n == 30, load hybrid_func2_M_D30,
%             elseif n == 50, load hybrid_func2_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 20
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func2_data
%             if n == 2, load hybrid_func2_M_D2,
%             elseif n == 10, load hybrid_func2_M_D10,
%             elseif n == 30, load hybrid_func2_M_D30,
%             elseif n == 50, load hybrid_func2_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 21
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func3_data
%             if n == 2, load hybrid_func3_M_D2,
%             elseif n == 10, load hybrid_func3_M_D10,
%             elseif n == 30, load hybrid_func3_M_D30,
%             elseif n == 50, load hybrid_func3_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 22
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func3_data
%             if n == 2, load hybrid_func3_HM_D2,
%             elseif n == 10, load hybrid_func3_HM_D10,
%             elseif n == 30, load hybrid_func3_HM_D30,
%             elseif n == 50, load hybrid_func3_HM_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 23
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func3_data
%             if n == 2, load hybrid_func3_M_D2,
%             elseif n == 10, load hybrid_func3_M_D10,
%             elseif n == 30, load hybrid_func3_M_D30,
%             elseif n == 50, load hybrid_func3_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 24
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func4_data
%             if n == 2, load hybrid_func4_M_D2,
%             elseif n == 10, load hybrid_func4_M_D10,
%             elseif n == 30, load hybrid_func4_M_D30,
%             elseif n == 50, load hybrid_func4_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 25
% 
%             lu = [2 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func4_data
%             if n == 2, load hybrid_func4_M_D2,
%             elseif n == 10, load hybrid_func4_M_D10,
%             elseif n == 30, load hybrid_func4_M_D30,
%             elseif n == 50, load hybrid_func4_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
%     end