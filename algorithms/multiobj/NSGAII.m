%**************************************************************************
% References
    % [1] *Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal, and T. Meyarivan*, |A Fast
    % Elitist Multiobjective Genetic Algorithm: NSGA-II|, IEEE Transactions on 
    % Evolutionary Computation 6 (2002), no. 2, 182 ~ 197.
    % [2] *N. Srinivas and Kalyanmoy Deb*, |Multiobjective Optimization Using 
    % Nondominated Sorting in Genetic Algorithms|, Evolutionary Computation 2 
    % (1994), no. 3, 221 ~ 248.
% Author:     unclear
% Adapter : Algori
% Email : zmdsn@126.com
% programmed: Sept 29, 2016
%**************************************************************************
function [best_solution,best_fvalue,run_series,run_info,user_set] = NSGAII(func_name , algo_config,func_config,user_config)

    format long;
    format compact;
    rand('seed', sum(100*clock));
    mytime = cputime;
    global objfunc;
    objfunc = feval(func_name);
    % **********************初始化算法参数********************************** 
    if ~exist('algo_config','var')
        algo_config = [];  % 参数初始化
    end
    % 种群规模
    if isfield(algo_config,'PopuSize')
        pop = algo_config.PopuSize;
    else
        pop = 100;
    end
    % 测试函数运行次数
    if isfield(algo_config,'maxFES')
        maxFES = algo_config.maxFES;
    else
        maxFES = 60000;
    end
    % 算法迭代代数
    if isfield(algo_config,'Gene')
        gen = algo_config.Gene;
    else
        gen = 500;
    end
    M = objfunc.Nobj;
    V = objfunc.Xdim;
    min_range = objfunc.Xmin;
    max_range = objfunc.Xmax;
    
    %% Initialize the population
    % Population is initialized with random values which are within the
    % specified range. Each chromosome consists of the decision variables. Also
    % the value of the objective functions, rank and crowding distance
    % information is also added to the chromosome vector but only the elements
    % of the vector which has the decision variables are operated upon to
    % perform the genetic operations like corssover and mutation.
    chromosome = initialize_variables(pop, M, V, min_range, max_range);
    %% Sort the initialized population
    % Sort the population using non-domination-sort. This returns two columns
    % for each individual which are the rank and the crowding distance
    % corresponding to their position in the front they belong. At this stage
    % the rank and the crowding distance for each chromosome is added to the
    % chromosome vector for easy of computation.
    chromosome = non_domination_sort_mod(chromosome, M, V);

    %% Start the evolution process
    % The following are performed in each generation
    % * Select the parents which are fit for reproduction
    % * Perfrom crossover and Mutation operator on the selected parents
    % * Perform Selection from the parents and the offsprings
    % * Replace the unfit individuals with the fit individuals to maintain a
    %   constant population size.
    gtime = cputime;
    for i = 1 : gen
        
        % Select the parents
        % Parents are selected for reproduction to generate offspring. The
        % original NSGA-II uses a binary tournament selection based on the
        % crowded-comparision operator. The arguments are 
        % pool - size of the mating pool. It is common to have this to be half the
        %        population size.
        % tour - Tournament size. Original NSGA-II uses a binary tournament
        %        selection, but to see the effect of tournament size this is kept
        %        arbitary, to be choosen by the user.
        pool = round(pop/2);
        tour = 2;
        % Selection process
        % A binary tournament selection is employed in NSGA-II. In a binary
        % tournament selection process two individuals are selected at random
        % and their fitness is compared. The individual with better fitness is
        % selcted as a parent. Tournament selection is carried out until the
        % pool size is filled. Basically a pool size is the number of parents
        % to be selected. The input arguments to the function
        % tournament_selection are chromosome, pool, tour. The function uses
        % only the information from last two elements in the chromosome vector.
        % The last element has the crowding distance information while the
        % penultimate element has the rank information. Selection is based on
        % rank and if individuals with same rank are encountered, crowding
        % distance is compared. A lower rank and higher crowding distance is
        % the selection criteria.
        parent_chromosome = tournament_selection(chromosome, pool, tour);
        % Perfrom crossover and Mutation operator
        % The original NSGA-II algorithm uses Simulated Binary Crossover (SBX) and
        % Polynomial  mutation. Crossover probability pc = 0.9 and mutation
        % probability is pm = 1/n, where n is the number of decision variables.
        % Both real-coded GA and binary-coded GA are implemented in the original
        % algorithm, while in this program only the real-coded GA is considered.
        % The distribution indeices for crossover and mutation operators as mu = 20
        % and mum = 20 respectively.
        mu = 20;
        mum = 20;
        offspring_chromosome = genetic_operator(parent_chromosome,M, V, mu, mum, min_range, max_range);
        % Intermediate population
        % Intermediate population is the combined population of parents and
        % offsprings of the current generation. The population size is two
        % times the initial population.
        [main_pop,temp] = size(chromosome);
        [offspring_pop,temp] = size(offspring_chromosome);
        % temp is a dummy variable.
        clear temp
        % intermediate_chromosome is a concatenation of current population and the offspring population.
        intermediate_chromosome(1:main_pop,:) = chromosome;
        intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = offspring_chromosome;

        % Non-domination-sort of intermediate population
        % The intermediate population is sorted again based on non-domination sort
        % before the replacement operator is performed on the intermediate
        % population.
        intermediate_chromosome = non_domination_sort_mod(intermediate_chromosome, M, V);
        % Perform Selection
        % Once the intermediate population is sorted only the best solution is
        % selected based on it rank and crowding distance. Each front is filled in
        % ascending order until the addition of population size is reached. The
        % last front is included in the population based on the individuals with
        % least crowding distance
        
        chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
        %if ~mod(i,100)
        %    clc
        %    fprintf('%d generations completed\n',i);
        %end
%         cputime - gtime;
    end
    cputime - gtime;
   %% Result  % For each chromosome perform the following
    [mm,nn] = size(chromosome);
    parsolution = chromosome(:,1:V);
    myres = zeros(mm,M);
    for i=1:mm
        myres(i,:) = evaluate_objective(parsolution(i,:), M, V);
    end
    % Save the result in ASCII text format.
%     save nsga2x.txt  parsolution -ASCII
%     save nsga2pf.txt myres -ASCII

    %% Visualize
%     if M == 2
%         plot(chromosome(:,V + 1),chromosome(:,V + 2),'*');
%     elseif M ==3
%         plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'*');
%     end
    objfunc.multi_draw(chromosome(:,V+1:V+M)',M);
%     run_series = outcome;
    
    % END: function outcome =HypeGA(Xmin, Xmax)
%     disp(['      最优值  =  ' num2str(globalBestValue,'%10.5e') ]);
    bRGAtime = cputime - mytime;
    run_info = bRGAtime;
    display(['      NSGAII time = ' num2str(bRGAtime)]);
end


%% evaluate_objective start
function f = evaluate_objective(x, M, V)
     global objfunc;
     f = objfunc.obj_value(x(1:V)')';
end % evaluate_objective end

%% initialize_variables start
function f = initialize_variables(N, M, V, min_range, max_range)
    %% function f = initialize_variables(N, M, V, min_tange, max_range) 
    % This function initializes the chromosomes. Each chromosome has the
    % following at this stage
    %       * set of decision variables
    %       * objective function values
    % where,
    % N - Population size
    % M - Number of objective functions
    % V - Number of decision variables
    % min_range - A vector of decimal values which indicate the minimum value for each decision variable.
    % max_range - Vector of maximum possible values for decision variables.
    min = min_range;     max = max_range;
    % K is the total number of array elements. For ease of computation decision
    % variables and objective functions are concatenated to form a single
    % array. For crossover and mutation only the decision variables are used
    % while for selection, only the objective variable are utilized.
    K = M + V;

    %% Initialize each chromosome
    % For each chromosome perform the following (N is the population size)
    for i = 1 : N
        % Initialize the decision variables based on the minimum and maximum
        % possible values. V is the number of decision variable. A random
        % number is picked between the minimum and maximum possible values for
        % the each decision variable.
        for j = 1 : V
            f(i,j) = min(j) + (max(j) - min(j))*rand(1);
        end
        % For ease of computation and handling data the chromosome also has the
        % vlaue of the objective function concatenated at the end. The elements
        % V + 1 to K has the objective function valued. 
        % The function evaluate_objective takes one chromosome at a time,
        % infact only the decision variables are passed to the function along
        % with information about the number of objective functions which are
        % processed and returns the value for the objective functions. These
        % values are now stored at the end of the chromosome itself.
        
        f(i,V + 1: K) = evaluate_objective(f(i,:), M, V);
    end
end

%% genetic_operator start
function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit)
    %% function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit) 
    % This function is utilized to produce offsprings from parent chromosomes.
    % The genetic operators corssover and mutation which are carried out with
    % slight modifications from the original design. For more information read the document enclosed. 
    % parent_chromosome - the set of selected chromosomes.
    % M - number of objective functions
    % V - number of decision varaiables
    % mu - distribution index for crossover (read the enlcosed pdf file)
    % mum - distribution index for mutation (read the enclosed pdf file)
    % l_limit - a vector of lower limit for the corresponding decsion variables
    % u_limit - a vector of upper limit for the corresponding decsion variables
    % The genetic operation is performed only on the decision variables, that
    % is the first V elements in the chromosome vector. 

    [N,m] = size(parent_chromosome);
    clear m
    p = 1;
    % Flags used to set if crossover and mutation were actually performed. 
    was_crossover = 0;
    was_mutation = 0;

    for i = 1 : N
        % With 90 % probability perform crossover
        if rand(1) < 0.9
            % Initialize the children to be null vector.
            child_1 = [];
            child_2 = [];
            % Select the first parent
            parent_1 = round(N*rand(1));
            if parent_1 < 1
                parent_1 = 1;
            end
            % Select the second parent
            parent_2 = round(N*rand(1));
            if parent_2 < 1
                parent_2 = 1;
            end
            % Make sure both the parents are not the same. 
            while isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:))
                parent_2 = round(N*rand(1));
                if parent_2 < 1
                    parent_2 = 1;
                end
            end
            % Get the chromosome information for each randomnly selected parents
            parent_1 = parent_chromosome(parent_1,:);
            parent_2 = parent_chromosome(parent_2,:);
            % Perform corssover for each decision variable in the chromosome.
            for j = 1 : V
                % SBX (Simulated Binary Crossover).
                % For more information about SBX refer the enclosed pdf file.
                % Generate a random number
                u(j) = rand(1);
                if u(j) <= 0.5
                    bq(j) = (2*u(j))^(1/(mu+1));
                else
                    bq(j) = (1/(2*(1 - u(j))))^(1/(mu+1));
                end
                % Generate the jth element of first child
                child_1(j) = ...
                    0.5*(((1 + bq(j))*parent_1(j)) + (1 - bq(j))*parent_2(j));
                % Generate the jth element of second child
                child_2(j) = ...
                    0.5*(((1 - bq(j))*parent_1(j)) + (1 + bq(j))*parent_2(j));
                % Make sure that the generated element is within the specified
                % decision space else set it to the appropriate extrema.
                if child_1(j) > u_limit(j)
                    child_1(j) = u_limit(j);
                elseif child_1(j) < l_limit(j)
                    child_1(j) = l_limit(j);
                end
                if child_2(j) > u_limit(j)
                    child_2(j) = u_limit(j);
                elseif child_2(j) < l_limit(j)
                    child_2(j) = l_limit(j);
                end
            end
            % Evaluate the objective function for the offsprings and as before
            % concatenate the offspring chromosome with objective value.
            child_1(:,V + 1: M + V) = evaluate_objective(child_1, M, V);
            child_2(:,V + 1: M + V) = evaluate_objective(child_2, M, V);
            % Set the crossover flag. When crossover is performed two children
            % are generated, while when mutation is performed only only child is
            % generated.
            was_crossover = 1;
            was_mutation = 0;
        % With 10 % probability perform mutation. Mutation is based on polynomial mutation. 
        else
            % Select at random the parent.
            parent_3 = round(N*rand(1));
            if parent_3 < 1
                parent_3 = 1;
            end
            % Get the chromosome information for the randomnly selected parent.
            child_3 = parent_chromosome(parent_3,:);
            % Perform mutation on eact element of the selected parent.
            for j = 1 : V
               r(j) = rand(1);
               if r(j) < 0.5
                   delta(j) = (2*r(j))^(1/(mum+1)) - 1;
               else
                   delta(j) = 1 - (2*(1 - r(j)))^(1/(mum+1));
               end
               % Generate the corresponding child element.
               child_3(j) = child_3(j) + delta(j);
               % Make sure that the generated element is within the decision
               % space.
               if child_3(j) > u_limit(j)
                   child_3(j) = u_limit(j);
               elseif child_3(j) < l_limit(j)
                   child_3(j) = l_limit(j);
               end
            end
            % Evaluate the objective function for the offspring and as before
            % concatenate the offspring chromosome with objective value.    
            child_3(:,V + 1: M + V) = evaluate_objective(child_3, M, V);
            % Set the mutation flag
            was_mutation = 1;
            was_crossover = 0;
        end
        % Keep proper count and appropriately fill the child variable with all
        % the generated children for the particular generation.
        if was_crossover
            child(p,:) = child_1;
            child(p+1,:) = child_2;
            was_cossover = 0;
            p = p + 2;
        elseif was_mutation
            child(p,:) = child_3(1,1 : M + V);
            was_mutation = 0;
            p = p + 1;
        end
    end
    f = child;
end

%% replace_chromosome start
function f  = replace_chromosome(intermediate_chromosome, M, V,pop)
%% function f  = replace_chromosome(intermediate_chromosome,pro,pop)
% This function replaces the chromosomes based on rank and crowding
% distance. Initially until the population size is reached each front is
% added one by one until addition of a complete front which results in
% exceeding the population size. At this point the chromosomes in that
% front is added subsequently to the population based on crowding distance.

    [N, m] = size(intermediate_chromosome);

    % Get the index for the population sort based on the rank
    [temp,index] = sort(intermediate_chromosome(:,M + V + 1));
    clear temp m

    % Now sort the individuals based on the index
    for i = 1 : N
        sorted_chromosome(i,:) = intermediate_chromosome(index(i),:);
    end

    % Find the maximum rank in the current population
    max_rank = max(intermediate_chromosome(:,M + V + 1));

    % Start adding each front based on rank and crowing distance until the
    % whole population is filled.
    previous_index = 0;
    for i = 1 : max_rank
        % Get the index for current rank i.e the last the last element in the
        % sorted_chromosome with rank i. 
        current_index = max(find(sorted_chromosome(:,M + V + 1) == i));
        % Check to see if the population is filled if all the individuals with
        % rank i is added to the population. 
        if current_index > pop
            % If so then find the number of individuals with in with current
            % rank i.
            remaining = pop - previous_index;
            % Get information about the individuals in the current rank i.
            temp_pop = ...
                sorted_chromosome(previous_index + 1 : current_index, :);
            % Sort the individuals with rank i in the descending order based on
            % the crowding distance.
            [temp_sort,temp_sort_index] = ...
                sort(temp_pop(:, M + V + 2),'descend');
            % Start filling individuals into the population in descending order
            % until the population is filled.
            for j = 1 : remaining
                f(previous_index + j,:) = temp_pop(temp_sort_index(j),:);
            end
            return;
        elseif current_index < pop
            % Add all the individuals with rank i into the population.
            f(previous_index + 1 : current_index, :) = ...
                sorted_chromosome(previous_index + 1 : current_index, :);
        else
            % Add all the individuals with rank i into the population.
            f(previous_index + 1 : current_index, :) = ...
                sorted_chromosome(previous_index + 1 : current_index, :);
            return;
        end
        % Get the index for the last added individual.
        previous_index = current_index;
    end
end

%% tournament_selection start
function f = tournament_selection(chromosome, pool_size, tour_size)
%% function tournament_selection(chromosome, pool_size, tour_size) 
% is the selection policy for selecting the individuals for the mating 
% pool. The selection is based on tournament selection. Argument 
% |chromosome| is the current generation population from which the 
% individuals are selected to form a mating pool of size |pool_size| after 
% performing tournament selection, with size of the tournament being 
% |tour_size|. By varying the tournament size the selection pressure can be
% adjusted. But for NSGA-II the tour_size is fixed to two, but the user may
% feel free to experiment with different tournament size. Also it has been
% observed that a tournament size of more than five has no significant
% meaning. 
%% Tournament selection process
% In a tournament selection process n individuals are selected at random,
% where n is equal to |tour_size|. From these individuals only one is selected
% and is added to the mating pool, where size of the mating pool is
% |pool_size|. Selection is performed based on two criteria. First and
% foremost is the rank or the front in which the solutions reside.
% Individuals with lower rank are selected. Secondly if the rank of two
% individuals are the same then, the crowding distance is compared.
% Individuals with greater crowding distance is selcted. 

    % Get the size of chromosome. The number of chromosome is not important
    % while the number of elements in chromosome are important.
    [pop, variables] = size(chromosome);
    % The peunltimate element contains the information about rank.
    rank = variables - 1;
    % The last element contains information about crowding distance.
    distance = variables;

    % Until the mating pool is filled, perform tournament selection
    for i = 1 : pool_size
        % Select n individuals at random, where n = tour_size
        for j = 1 : tour_size
            % Select an individual at random
            candidate(j) = round(pop*rand(1));
            % Make sure that the array starts from one. 
            if candidate(j) == 0
                candidate(j) = 1;
            end
            if j > 1
                % Make sure that same candidate is not choosen.
                while ~isempty(find(candidate(1 : j - 1) == candidate(j)))
                    candidate(j) = round(pop*rand(1));
                    if  candidate(j) == 0
                        candidate(j) = 1;
                    end
                end
            end
        end
        % Collect information about the selected candidates.
        for j = 1 : tour_size
            c_obj_rank(j) = chromosome(candidate(j),rank);
            c_obj_distance(j) = chromosome(candidate(j),distance);
        end
        % Find the candidate with the least rank
        min_candidate = find(c_obj_rank == min(c_obj_rank));
        % If more than one candiate have the least rank then find the candidate
        % within that group having the maximum crowding distance.
        if length(min_candidate) ~= 1
            max_candidate = ...
            find(c_obj_distance(min_candidate) == max(c_obj_distance(min_candidate)));
            % If a few individuals have the least rank and have maximum crowding
            % distance, select only one individual (not at random). 
            if length(max_candidate) ~= 1
                max_candidate = max_candidate(1);
            end
            % Add the selected individual to the mating pool
            f(i,:) = chromosome(candidate(min_candidate(max_candidate)),:);
        else
            % Add the selected individual to the mating pool
            f(i,:) = chromosome(candidate(min_candidate(1)),:);
        end
    end
end

function f = non_domination_sort_mod(x, M, V)
    %% function f = non_domination_sort_mod(x, M, V)
    % This function sort the current popultion based on non-domination. All the
    % individuals in the first front are given a rank of 1, the second front
    % individuals are assigned rank 2 and so on. After assigning the rank the
    % crowding in each front is calculated.
    [N, m] = size(x);   clear m
    % Initialize the front number to 1.
    front = 1; 
    % There is nothing to this assignment, used only to manipulate easily in MATLAB.
    F(front).f = [];   individual = [];
    %% Non-Dominated sort. 
    % The initialized population is sorted based on non-domination. The fast
    % sort algorithm [1] is described as below for each
    % ?for each individual p in main population P do the following
    %   ?Initialize Sp = []. This set would contain all the individuals that is
    %     being dominated by p.
    %   ?Initialize np = 0. This would be the number of individuals that dominate p.
    %   ?for each individual q in P
    %       * if p dominated q then
    %           ?add q to the set Sp i.e. Sp = Sp ? {q}
    %       * else if q dominates p then
    %           ?increment the domination counter for p i.e. np = np + 1
    %   ?if np = 0 i.e. no individuals dominate p then p belongs to the first
    %     front; Set rank of individual p to one i.e prank = 1. Update the first
    %     front set by adding p to front one i.e F1 = F1 ? {p}
    % ?This is carried out for all the individuals in main population P.
    % ?Initialize the front counter to one. i = 1
    % ?following is carried out while the ith front is nonempty i.e. Fi != []
    %   ?Q = []. The set for storing the individuals for (i + 1)th front.
    %   ?for each individual p in front Fi
    %       * for each individual q in Sp (Sp is the set of individuals
    %         dominated by p)
    %           ?nq = nq?1, decrement the domination count for individual q.
    %           ?if nq = 0 then none of the individuals in the subsequent
    %             fronts would dominate q. Hence set qrank = i + 1. Update
    %             the set Q with individual q i.e. Q = Q ? q.
    %   ?Increment the front counter by one.
    %   ?Now the set Q is the next front and hence Fi = Q.
    %
    % This algorithm is better than the original NSGA ([2]) since it utilize
    % the informatoion about the set that an individual dominate (Sp) and
    % number of individuals that dominate the individual (np).
    for i = 1 : N
        % Number of individuals that dominate this individual
        individual(i).n = 0;
        % Individuals which this individual dominate
        individual(i).p = [];
        for j = 1 : N
            dom_less = 0;
            dom_equal = 0;
            dom_more = 0;
            for k = 1 : M
                if (x(i,V + k) < x(j,V + k))
                    dom_less = dom_less + 1;
                elseif (x(i,V + k) == x(j,V + k))
                    dom_equal = dom_equal + 1;
                else
                    dom_more = dom_more + 1;
                end
            end
            if dom_less == 0 && dom_equal ~= M
                individual(i).n = individual(i).n + 1;
            elseif dom_more == 0 && dom_equal ~= M
                individual(i).p = [individual(i).p j];
            end
        end   
        if individual(i).n == 0
            x(i,M + V + 1) = 1;
            F(front).f = [F(front).f i];
        end
    end
    % Find the subsequent fronts
    while ~isempty(F(front).f)
       Q = [];
       for i = 1 : length(F(front).f)
           if ~isempty(individual(F(front).f(i)).p)
                for j = 1 : length(individual(F(front).f(i)).p)
                    individual(individual(F(front).f(i)).p(j)).n = ...
                        individual(individual(F(front).f(i)).p(j)).n - 1;
                    if individual(individual(F(front).f(i)).p(j)).n == 0
                        x(individual(F(front).f(i)).p(j),M + V + 1) = ...
                            front + 1;
                        Q = [Q individual(F(front).f(i)).p(j)];
                    end
                end
           end
       end
       front =  front + 1;
       F(front).f = Q;
    end
    [temp,index_of_fronts] = sort(x(:,M + V + 1));
    for i = 1 : length(index_of_fronts)
        sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
    end
    current_index = 0;
    %% Crowding distance
    %The crowing distance is calculated as below
    % ?For each front Fi, n is the number of individuals.
    %   ?initialize the distance to be zero for all the individuals i.e. Fi(dj ) = 0,
    %     where j corresponds to the jth individual in front Fi.
    %   ?for each objective function m
    %       * Sort the individuals in front Fi based on objective m i.e. I =
    %         sort(Fi,m).
    %       * Assign infinite distance to boundary values for each individual
    %         in Fi i.e. I(d1) = ? and I(dn) = ?
    %       * for k = 2 to (n ? 1)
    %           ?I(dk) = I(dk) + (I(k + 1).m ? I(k ? 1).m)/fmax(m) - fmin(m)
    %           ?I(k).m is the value of the mth objective function of the kth
    %             individual in I

    % Find the crowding distance for each individual in each front
    for front = 1 : (length(F) - 1)
    %    objective = [];
        distance = 0;
        y = [];
        previous_index = current_index + 1;
        for i = 1 : length(F(front).f)
            y(i,:) = sorted_based_on_front(current_index + i,:);
        end
        current_index = current_index + i;
        % Sort each individual based on the objective
        sorted_based_on_objective = [];
        for i = 1 : M
            [sorted_based_on_objective, index_of_objectives] = ...
                sort(y(:,V + i));
            sorted_based_on_objective = [];
            for j = 1 : length(index_of_objectives)
                sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
            end
            f_max = ...
                sorted_based_on_objective(length(index_of_objectives), V + i);
            f_min = sorted_based_on_objective(1, V + i);
            y(index_of_objectives(length(index_of_objectives)),M + V + 1 + i)...
                = Inf;
            y(index_of_objectives(1),M + V + 1 + i) = Inf;
             for j = 2 : length(index_of_objectives) - 1
                next_obj  = sorted_based_on_objective(j + 1,V + i);
                previous_obj  = sorted_based_on_objective(j - 1,V + i);
                if (f_max - f_min == 0)
                    y(index_of_objectives(j),M + V + 1 + i) = Inf;
                else
                    y(index_of_objectives(j),M + V + 1 + i) = ...
                         (next_obj - previous_obj)/(f_max - f_min);
                end
             end
        end
        distance = [];
        distance(:,1) = zeros(length(F(front).f),1);
        for i = 1 : M
            distance(:,1) = distance(:,1) + y(:,M + V + 1 + i);
        end
        y(:,M + V + 2) = distance;
        y = y(:,1 : M + V + 2);
        z(previous_index:current_index,:) = y;
    end
    f = z();
    %% References
    % [1] *Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal, and T. Meyarivan*, |A Fast
    % Elitist Multiobjective Genetic Algorithm: NSGA-II|, IEEE Transactions on 
    % Evolutionary Computation 6 (2002), no. 2, 182 ~ 197.
    % [2] *N. Srinivas and Kalyanmoy Deb*, |Multiobjective Optimization Using 
    % Nondominated Sorting in Genetic Algorithms|, Evolutionary Computation 2 
    % (1994), no. 3, 221 ~ 248.
end

