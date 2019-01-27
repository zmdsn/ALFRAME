function [best_fvalue,best_solution,run_series,run_info,user_set] = SaDE(func_name , algo_config,func_config,user_config)

%****************************************************************************
%Reference:  
    % A. K. Qin, V. L. Huang, and P. N. Suganthan,°∞Differential evolution
    % algorithm with strategy adaptation for global numerical optimization,°±
    % IEEE Trans. Evolut. Comput., vol. 13, no. 2, pp. 398®C417, Apr. 2009.
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
    [algo_config, func_config, objfunc, NP, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
    %% ****************==- main body  -==***********************

    lu = [Xmin * ones(1, DIM); Xmax * ones(1, DIM)];     
    outcome = [];  % record the best results

    %% *****************==- Main body which was provided by the authors -==**********************
        numst = 4;
        aaaa = cell(1, numst);
        learngen = 50;

        lpcount = [];
        npcount = [];

        %Record the number of success or failure
        ns = [];
        nf = [];
        %Record the success rate
        pfit = ones(1, numst);
        %Record the median of CR
        ccm = 0.5 * ones(1, numst);

        %-----Initialize population and some arrays-------------------------------
        pop = zeros(NP, DIM); %initialize pop to gain speed
        XRRmin = repmat(lu(1, :), NP, 1);
        XRRmax = repmat(lu(2, :), NP, 1);
        pop = XRRmin + (XRRmax - XRRmin).* rand(NP, DIM);

        popold   = zeros(size(pop));   % toggle population
        val    = zeros(1, NP);      % create and reset the "cost array"
        DE_gbest  = zeros(1, DIM);      % best population member ever
        nfeval = 0;           % number of function evaluations

        %------Evaluate the best member after initialization----------------------

        val  = f(pop',func_config);
        [DE_gbestval, ibest] = min(val);
        DE_gbest = pop(ibest, :);

        %     ibest  = 1;            % start with first population member
        %     val(1)  = benchmark_func(pop(ibest, :), problem, o, A, M, a, alpha, b);
        %     DE_gbestval = val(1);         % best objective function value so far
        %     nfeval  = nfeval + 1;
        %     for i = 2:NP            % check the remaining members
        %       val(i) = benchmark_func(pop(i, :), problem, o, A, M, a, alpha, b);
        %       nfeval  = nfeval + 1;
        %       if (val(i) < DE_gbestval)      % if member is better
        %         ibest  = i;         % save its location
        %         DE_gbestval = val(i);
        %       end
        %     end
        %     DE_gbest = pop(ibest, :);     % best member of current iteration

        %------DE-Minimization---------------------------------------------
        %------popold is the population which has to compete. It is--------
        %------static through one iteration. pop is the newly--------------
        %------emerging population.----------------------------------------

        pm1 = zeros(NP, DIM);        % initialize population matrix 1
        pm2 = zeros(NP, DIM);        % initialize population matrix 2
        pm3 = zeros(NP, DIM);        % initialize population matrix 3
        pm4 = zeros(NP, DIM);        % initialize population matrix 4
        pm5 = zeros(NP, DIM);        % initialize population matrix 5
        bm  = zeros(NP, DIM);        % initialize DE_gbestber matrix
        ui  = zeros(NP, DIM);        % intermediate population of perturbed vectors
        mui = zeros(NP, DIM);        % mask for intermediate population
        mpo = zeros(NP, DIM);        % mask for old population
        rot = (0 : 1 : NP-1);        % rotating index array (size NP)
        rt  = zeros(NP);         % another rotating index array
        a1  = zeros(NP);         % index array
        a2  = zeros(NP);         % index array
        a3  = zeros(NP);         % index array
        a4  = zeros(NP);         % index array
        a5  = zeros(NP);         % index array
        ind = zeros(4);

        iter = 1;
        nfeval = NP;
        while nfeval < maxFES

            popold = pop;             % save the old population
            ind = randperm(4);        % index pointer array
            a1  = randperm(NP);       % shuffle locations of vectors
            rt = rem(rot + ind(1), NP);     % rotate indices by ind(1) positions
            a2  = a1(rt + 1);         % rotate vector locations
            rt = rem(rot + ind(2), NP);
            a3  = a2(rt + 1);
            rt = rem(rot + ind(3), NP);
            a4  = a3(rt + 1);
            rt = rem(rot + ind(4), NP);
            a5  = a4(rt + 1);

            pm1 = popold(a1, :);       % shuffled population 1
            pm2 = popold(a2, :);       % shuffled population 2
            pm3 = popold(a3, :);       % shuffled population 3
            pm4 = popold(a4, :);       % shuffled population 4
            pm5 = popold(a5, :);       % shuffled population 5

            bm = repmat(DE_gbest, NP, 1);

            if (iter >= learngen)
                for i = 1:numst
                    if  ~isempty(aaaa{i})
                        ccm(i) = median(aaaa{i}(:, 1));
                        d_index = find(aaaa{i}(:, 2) == aaaa{i}(1, 2));
                        aaaa{i}(d_index, :) = [];
                    else
                        ccm(i) = rand;
                    end
                end
            end

            for i = 1 : numst
                cc_tmp = [];
                for k = 1 : NP
                    tt = normrnd(ccm(i), 0.1);
                    while tt > 1 | tt < 0
                        tt = normrnd(ccm(i), 0.1);
                    end
                    cc_tmp = [cc_tmp; tt];
                end
                cc(:, i) = cc_tmp;
            end

            % Stochastic universal sampling
            rr = rand;
            spacing = 1/NP;
            randnums = sort(mod(rr : spacing : 1 + rr - 0.5 * spacing, 1));

            normfit = pfit / sum(pfit);
            partsum = 0;
            count(1) = 0;
            stpool = [];

            for i = 1 : length(pfit)
                partsum = partsum + normfit(i);
                count(i + 1) = length(find(randnums < partsum));
                select(i, 1) = count(i + 1) - count(i);
                stpool = [stpool; ones(select(i, 1), 1) * i];
            end
            stpool = stpool(randperm(NP));

            for i = 1 : numst
                atemp = zeros(1, NP);
                aaa{i} = atemp;
                index{i} = [];
                if ~isempty(find(stpool == i))
                    index{i} = find(stpool == i);
                    atemp(index{i}) = 1;
                    aaa{i} = atemp;
                end
            end

            aa = zeros(NP, DIM);
            for i = 1 : numst
                 % all random numbers < CR are 1, 0 otherwise
                aa(index{i}, :) = rand(length(index{i}), DIM) < repmat(cc(index{i}, i), 1, DIM);
            end
            mui = aa;

            % jrand
            dd = ceil(DIM * rand(NP, 1));
            for kk = 1 : NP
                mui(kk, dd(kk)) = 1;
            end
            mpo = mui < 0.5;         % inverse mask to mui

            for i = 1 : numst
                %-----------jitter---------
                F = [];
                m = length(index{i});
                F = normrnd(0.5, 0.3, m, 1);
                F = repmat(F, 1, DIM);
                if i == 1
                    ui(index{i}, :) = pm3(index{i}, :) + F .* (pm1(index{i}, :) - pm2(index{i}, :));     % differential variation
                    ui(index{i}, :) = popold(index{i}, :) .* mpo(index{i}, :) + ui(index{i}, :) .* mui(index{i}, :);   % crossover
                end
                if i == 2
                    ui(index{i}, :) = popold(index{i}, :) + F .* (bm(index{i}, :)-popold(index{i}, :)) + F .* (pm1(index{i}, :) - pm2(index{i}, :) + pm3(index{i}, :) - pm4(index{i}, :));    % differential variation
                    ui(index{i}, :) = popold(index{i}, :) .* mpo(index{i}, :) + ui(index{i}, :) .* mui(index{i}, :);   % crossover
                end
                if i == 3
                    ui(index{i}, :) = pm5(index{i}, :) + F .* (pm1(index{i}, :) - pm2(index{i}, :) + pm3(index{i}, :) - pm4(index{i}, :));    % differential variation
                    ui(index{i}, :) = popold(index{i}, :) .* mpo(index{i}, :) + ui(index{i}, :) .* mui(index{i}, :);   % crossover
                end
                if i == 4
                    ui(index{i}, :) = popold(index{i}, :) + rand .* (pm5(index{i}, :)-popold(index{i}, :)) + F .* (pm1(index{i}, :) - pm2(index{i}, :));
                end
            end

            for i = 1 : NP
                outbind = find(ui(i, :) < lu(1, :));
                XRmin = lu(1, :);
                XRmax = lu(2, :);
                if size(outbind, 2) ~= 0
                    ui(i, outbind) = XRmin(outbind) + (XRmax(outbind) - XRmin(outbind)) .* rand(1, size(outbind, 2));
                end
                outbind = find(ui(i, :) > lu(2, :));
                if size(outbind, 2) ~= 0
                    ui(i, outbind) = XRmin(outbind) + (XRmax(outbind) - XRmin(outbind)) .* rand(1, size(outbind, 2));
                end
            end

            lpcount = zeros(1, numst);
            npcount = zeros(1, numst);

            tempval = f(ui',func_config);  % check cost of competitor
            nfeval  = nfeval + NP;

            for i = 1 : NP

                if (tempval(i) <= val(i)) % if competitor is better than value in "cost array"

                    pop(i, :) = ui(i, :);  % replace old vector with new one (for new iteration)
                    val(i)  = tempval(i);  % save value in "cost array"

                    tlpcount = zeros(1, numst);
                    for j = 1 : numst
                        temp = aaa{j};
                        tlpcount(j) = temp(i);
                        if tlpcount(j) == 1
                            aaaa{j} = [aaaa{j}; cc(i, j) iter];
                        end
                    end
                    lpcount = [lpcount; tlpcount];

                else

                    tnpcount = zeros(1, numst);
                    for j = 1:numst
                        temp = aaa{j};
                        tnpcount(j) = temp(i);
                    end
                    npcount = [npcount; tnpcount];

                end

            end %---end for imember = 1:NP

            ns = [ns; sum(lpcount, 1)];
            nf = [nf; sum(npcount, 1)];

            if iter >= learngen,
                for i = 1 : numst
                    if (sum(ns(:, i)) + sum(nf(:, i))) == 0
                        pfit(i) = 0.01;
                    else
                        pfit(i) = sum(ns(:, i)) / (sum(ns(:, i)) + sum(nf(:, i))) + 0.01;
                    end
                end
                if ~isempty(ns), ns(1, :) = [];  end
                if ~isempty(nf), nf(1, :) = [];  end
            end
            iter = iter + 1;
            outcome = [outcome; [nfeval min(val)]];
        end
    %% ****************==- collating the results -==*********************
    [best_fvalue,minidx]= min(val);
    best_solution = pop(minidx, :)';
    run_series = outcome;
    Altime = cputime - mytime;                  % º∆À„ ±º‰
    run_info = [Altime];
end



