function [best_fvalue,best_solution,run_series,run_info,user_set] = jDE(func_name , algo_config,func_config,user_config)

%**************************************************************************************************
%Reference:  1) J. Brest, S. Greiner, B. Boskovic, M. Mernik, and V. Zumer, ��Self-adapting
%                     control parameters in differential evolution: A comparative study on numerical
%                     benchmark problems,�� IEEE Trans. Evolut. Comput., vol. 10, no. 6,
%                     pp. 646�C657, Dec. 2006.
%                     2) J. Zhang and A. C. Sanderson, ��JADE: adaptive differential evolution
%                     with optional external archive,�� IEEE Trans. Evolut. Comput., vol. 13,
%                     no. 5, pp. 945-958, 2009.
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
    default_set = struct('PopuSize' , 100,...
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
     
    tau1 = 0.1;
    tau2 = 0.1;

    F = 0.5 * ones(popsize, 1);
    CR = 0.9 * ones(popsize, 1);
    lu = [Xmin * ones(1, DIM); Xmax * ones(1, DIM)];    
    outcome = [];  % record the best results
    
    %Main body which was provided by Dr. J. Zhang
%     time = 1;
%     % The total number of runs
%     totalTime = 1;
%     while time <= totalTime
        
        % Initialize the main population
        popold = repmat(lu(1, :), popsize, 1) + rand(popsize, DIM) .* (repmat(lu(2, :)-lu(1, :), popsize, 1));
        valParents =  f(popold',func_config)';
        FES = popsize;
        outcome = [outcome; [FES min(valParents)]];
        while FES < maxFES
            pop = popold;      % the old population becomes the current population
            Fold = F;
            CRold = CR;

            IF  = rand(popsize, 1) < tau1;
            ICR = rand(popsize, 1) < tau2;

            F(IF)  = 0.1 + 0.9 * rand(sum(IF), 1);
            CR(ICR) = 0.0 + 1.0 * rand(sum(ICR), 1);

            r0 = [1:popsize];
            [r1, r2, r3] = gnR1R2R3(popsize, r0);

            %  == == == == == = Mutation == == == == == == == == =
            vi  = pop(r1, :) + F(:, ones(1, DIM)) .* (pop(r2, :) - pop(r3, :));

            %�˴�ԭʼ�����vi = boundConstraint(vi,lu);���ǲ����boundConstraint������÷�������JADE�Ҹĳ�����
            vi = boundConstraint(vi,pop, lu); 

            % == == == == =  Crossover == == == == =
             % mask is used to indicate which elements of ui comes from the parent
            mask = rand(popsize, DIM) > CR(:, ones(1, DIM)); 
            % choose one position where the element of ui doesn't come from the parent
            rows = (1:popsize)';
            cols = floor(rand(popsize, 1) * DIM) + 1;
            jrand = sub2ind([popsize DIM], rows, cols);
            mask(jrand) = false;
            ui = vi;
            ui(mask) = pop(mask);

            valOffspring =  f(ui',func_config)';
            FES = FES + popsize;
            
%             objfunc.draw_running(pop');
            %  == == == == == = Selection == == == == == =
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents, I] = min([valParents, valOffspring], [], 2);
            popold = pop;
            popold(I == 2, :) = ui(I == 2, :);
            F(I == 1) = Fold(I == 1);
            CR(I == 1) = CRold(I == 1);
            outcome = [outcome; [FES min(valParents)]];
        end
    %% ****************==- collating the results -==*********************
        best_fvalue = min(valParents);
        best_solution = pop(I(1),:)';
        run_series = outcome;
        Altime = cputime - mytime;                  % ����ʱ��
        run_info = [Altime];
        
        
%         display(['      jDE   ' num2str(cputime - mytime)]);
%         display(['      ����ֵ��'   num2str(min(valParents))]);
end

function [r1,r2,r3] = gnR1R2R3(NP1, r0)

    % gnA1A2 generate two column vectors r1 and r2 of size NP1 & NP2, respectively
    %       r1's elements are choosen from {1, 2, ..., NP1} & r1(i) ~= r0(i)
    %       r2's elements are choosen from {1, 2, ..., NP2} & r2(i) ~= r1(i) & r2(i) ~= r0(i)
    %
    % Call:
    %        [r1 r2 ...] = gnA1A2(NP1)      % r0 is set to be (1:NP1)'
    %        [r1 r2 ...] = gnA1A2(NP1, r0)  % r0 should be of length NP1
    %
    % Version: 2.1   Date: 2008/07/01
    % Written by Jingqiao Zhang (jingqiao@gmail.com)
    warning off all
    NP0 = length(r0);

    r1  = floor(rand(1, NP0) * NP1) + 1;
    for i = 1 : inf
        pos = (r1 == r0);
        if sum(pos) == 0
            break;
        else   % regenerate r1 if it is equal to r0
            r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r1 in 1000 iterations');
        end
    end

    r2  = floor(rand(1, NP0) * NP1) + 1;
    for i = 1 : inf
        pos = ((r2 == r1) | (r2 == r0));
        if sum(pos) == 0
            break;
        else   % regenerate r2 if it is equal to r0 or r1
            r2(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r2 in 1000 iterations');
        end
    end

    r3  = floor(rand(1, NP0) * NP1) + 1;
    for i = 1 : inf
        pos = ((r3 == r1) | (r3 == r0) | (r3 == r2));
        if sum(pos) == 0
            break;
        else   % regenerate r2 if it is equal to r0 or r1
            r3(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r3 in 1000 iterations');
        end
    end
 
end

function vi = boundConstraint (vi, pop, lu)

    % if the boundary constraint is violated, set the value to be the middle
    % of the previous value and the bound
    %
    % Version: 1.1   Date: 11/20/2007
    % Written by Jingqiao Zhang, jingqiao@gmail.com

    [NP, D] = size(pop);  % the population size and the problem's dimension

    %% check the lower bound
    xl = repmat(lu(1, :), NP, 1);
    pos = vi < xl;
    vi(pos) = (pop(pos) + xl(pos)) / 2;

    %% check the upper bound
    xu = repmat(lu(2, :), NP, 1);
    pos = vi > xu;
    vi(pos) = (pop(pos) + xu(pos)) / 2;
end
