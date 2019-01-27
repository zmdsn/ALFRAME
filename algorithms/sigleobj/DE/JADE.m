function [best_fvalue,best_solution,run_series,run_info,user_set] = JADE(func_name , algo_config,func_config,user_config)
%**************************************************************************************************
%Reference:  J. Zhang and A. C. Sanderson, "JADE: adaptive differential evolution
%                     with optional external archive," IEEE Trans. Evolut. Comput., vol. 13,
%                     no. 5, pp. 945-958, 2009.
%
% JADE with Archive
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
                         'maxFES'   , 600000  );
                     
    if ~exist('func_name','var') func_name = @sin; end
    if ~exist('algo_config','var') algo_config = []; end
    if ~exist('func_config','var') func_config = []; end
    % get parameters
    [algo_config, func_config, objfunc, popsize, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
    %% ****************==- main body  -==***********************
    n = DIM ;         % 获取目标函数类的一些参数
    lu = [objfunc.Xmin'; objfunc.Xmax'];

    outcome = [];
    time = 1;
    % Initialize the main population
    popold = repmat(lu(1, :), popsize, 1) + rand(popsize, n) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));

    valParents = f( popold' , func_config)'; % benchmark_func(, problem, o, A, M, a, alpha, b);
    c = 1/10;
    p = 0.05;
    % CR=0.95*ones(100,1);
    %   F=0.75*ones(100,1);

      

        CRm = 0.5;
        Fm = 0.5;
        CRsigma = 0.1;
        Fsigma = 0.1;
        Afactor = 1;

        archive.NP = Afactor * popsize; % the maximum size of the archive
        archive.pop = zeros(0, n); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions

        %% the values and indices of the best solutions
        [valBest, indBest] = sort(valParents, 'ascend');
        
        SF = [];
        SCR = [];        

        FES = 0;
        while FES < maxFES %& min(fit)>error_value(problem)

            [Fi,CRi] = randFCR(popsize, CRm, CRsigma, Fm,  Fsigma);
            pop = popold; % the old population becomes the current population
            r0 = [1 : popsize];
            popAll = [pop; archive.pop];
            [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);

            % Find the p-best solutions
            pNP = max(round(p * popsize), 2); % choose at least two best solutions
            randindex = ceil(rand(1, popsize) * pNP); % select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop(indBest(randindex), :);  % randomly choose one of the top 100p% solutions
            % == == == == == == == == == == == == == == == Mutation == == == == == == == == == == == == == 
            vi = pop + Fi(:, ones(1, n)) .* (pbest - pop) + Fi(:, ones(1, n)).* (pop(r1, :) - popAll(r2, :));   
            vi = boundConstraint(vi, pop, lu);
            % == == == == = Crossover == == == == =
            %if FES>100000
            mask = rand(popsize, n) > CRi(:, ones(1, n)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize)'; cols = floor(rand(popsize, 1) * n)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize n], rows, cols);
            mask(jrand) = false;
            ui = vi;
            ui(mask) = pop(mask);
            valOffspring =f( ui',func_config )'; % benchmark_func(ui, problem, o, A, M, a, alpha, b);
            FES = FES + popsize;
            
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents, I] = min([valParents, valOffspring], [], 2);
            popold = pop;
            archive = updateArchive(archive, popold(I == 2, :), valParents(I == 2));
            popold(I == 2, :) = ui(I == 2, :);
            SF = Fi(I == 2);
            SCR = CRi(I == 2);
            CRm = (1-c)*CRm + c*mean(SCR);
            Fm = (1-c)*Fm + c*meanL(SF);
            [valBest indBest] = sort(valParents, 'ascend');
            outcome = [outcome ;FES  valBest(1)];
        end

        
%         time = time + 1;
%     end

    %% ****************==- collating the results -==*********************
    best_fvalue = valBest(1);
    best_solution = pop(indBest(1), :)';
    run_series = outcome;
    Altime = cputime - mytime;                  % 计算时间
    run_info = [Altime];    
end


function [F,CR] = randFCR(NP, CRm, CRsigma, Fm,  Fsigma)

% this function generate CR according to a normal distribution with mean "CRm" and sigma "CRsigma"
%           If CR > 1, set CR = 1. If CR < 0, set CR = 0.
% this function generate F  according to a cauchy distribution with location parameter "Fm" and scale parameter "Fsigma"
%           If F > 1, set F = 1. If F <= 0, regenrate F.
%
% Version: 1.1   Date: 11/20/2007
% Written by Jingqiao Zhang (jingqiao@gmail.com)

    %% generate CR
    CR = CRm + CRsigma * randn(NP, 1);
    CR = min(1, max(0, CR));                % truncated to [0 1]

    %% generate F
    F = randCauchy(NP, 1, Fm, Fsigma);
    F = min(1, F);                          % truncation

    % we don't want F = 0. So, if F<=0, we regenerate F (instead of trucating it to 0)
    pos = find(F <= 0);
    while ~ isempty(pos)
        F(pos) = randCauchy(length(pos), 1, Fm, Fsigma);
        F = min(1, F);                      % truncation
        pos = find(F <= 0);
    end

end

% Cauchy distribution: cauchypdf = @(x, mu, delta) 1/pi*delta./((x-mu).^2+delta^2)
function result = randCauchy(m, n, mu, delta)

% http://en.wikipedia.org/wiki/Cauchy_distribution
% size(mu)
% size(delta * tan(pi * (rand(m, n) - 0.5)))
result = mu + delta * tan(pi * (rand(m, n) - 0.5));
end

function archive = updateArchive(archive, pop, funvalue)
% Update the archive with input solutions
%   Step 1: Add new solution to the archive
%   Step 2: Remove duplicate elements
%   Step 3: If necessary, randomly remove some solutions to maintain the archive size
%
% Version: 1.1   Date: 2008/04/02
% Written by Jingqiao Zhang (jingqiao@gmail.com)

    if archive.NP == 0, return; end

    if size(pop, 1) ~= size(funvalue,1), error('check it'); end

    % Method 2: Remove duplicate elements
    popAll = [archive.pop; pop ];
    funvalues = [archive.funvalues; funvalue ];
    [dummy IX]= unique(popAll, 'rows');
    if length(IX) < size(popAll, 1) % There exist some duplicate solutions
      popAll = popAll(IX, :);
      funvalues = funvalues(IX, :);
    end

    if size(popAll, 1) <= archive.NP   % add all new individuals
      archive.pop = popAll;
      archive.funvalues = funvalues;
    else                % randomly remove some solutions
      rndpos = randperm(size(popAll, 1)); % equivelent to "randperm";
      rndpos = rndpos(1 : archive.NP);

      archive.pop = popAll  (rndpos, :);
      archive.funvalues = funvalues(rndpos, :);
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

function [r1, r2] = gnR1R2(NP1, NP2, r0)
 warning('off');
% gnA1A2 generate two column vectors r1 and r2 of size NP1 & NP2, respectively
%    r1's elements are choosen from {1, 2, ..., NP1} & r1(i) ~= r0(i)
%    r2's elements are choosen from {1, 2, ..., NP2} & r2(i) ~= r1(i) & r2(i) ~= r0(i)
%
% Call:
%    [r1 r2 ...] = gnA1A2(NP1)   % r0 is set to be (1:NP1)'
%    [r1 r2 ...] = gnA1A2(NP1, r0) % r0 should be of length NP1
%
% Version: 2.1  Date: 2008/07/01
% Written by Jingqiao Zhang (jingqiao@gmail.com)

    NP0 = length(r0);

    r1 = floor(rand(1, NP0) * NP1) + 1;
    for i = 1 : inf
        pos = (r1 == r0);
        if sum(pos) == 0
            break;
        else % regenerate r1 if it is equal to r0
            r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r1 in 1000 iterations');
        end
    end

    r2 = floor(rand(1, NP0) * NP2) + 1;
    for i = 1 : inf
        pos = ((r2 == r1) | (r2 == r0));
        if sum(pos)==0
            break;
        else % regenerate r2 if it is equal to r0 or r1
            r2(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r2 in 1000 iterations');
        end
    end
end

% mean_L = \frac{\sum_{F \in S_F}F^2}{\sum_{F \in S_F}F}
function result = meanL(F)
    result = F'*F/sum(F);
end
