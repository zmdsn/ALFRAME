function [best_fvalue,best_solution,run_series,run_info,user_set] = FWA(func_name , algo_config,func_config,user_config)
    mytime = cputime;
    global total_fiteval_times;
    global Coef_Explosion_Amplitude;
    global Max_Sparks_Num;
    global Min_Sparks_Num;
    global Coef_Spark_Num;

    total_fiteval_times = 0;
    Coef_Explosion_Amplitude = 40;
    Max_Sparks_Num = 40;
    Min_Sparks_Num = 2;
    Coef_Spark_Num = 50;
    global objfunc ;
    %% *****************==- ��ʼ�����Ժ��� -==**********************
    if exist('func_config','var')
        objfunc = feval(func_name,func_config);  % ʵ��һ��Ŀ�꺯����
    else
        objfunc = feval(func_name);  % ʵ��һ��Ŀ�꺯����
        func_config = [];
    end
    objfunc = feval(func_name,func_config);
    
    %% *****************==- ��ʼ���㷨���� -==**********************
    if ~exist('algo_config','var')
        algo_config = [];  % �����ʼ��?
    end
    % ��Ⱥ��ģ
    if isfield(algo_config,'PopuSize')
        params.seednum = algo_config.PopuSize;
    else
        params.seednum = 10;
    end
    % ������
    if isfield(algo_config,'maxFES')
        params.maxEva = algo_config.maxFES;
    else
        params.maxEva = 60000;
    end

    outcome = [];
    params.dim              = objfunc.Xdim;        % problem dimension
    % params.seednum          = 5;
    params.sonnum           = 50; 
    % params.maxEva		    = 300000;        % the max iteration of the algorithm
    params.modStep          = 100;
    params.maxEva_mod100    = params.maxEva/params.modStep;
    params.gaussianNum		= 5;
    params.lowerBound =  objfunc.Xmin(1);
    params.upperBound = objfunc.Xmax(1);
    params.lowerInit =  objfunc.Xmin(1);
    params.upperInit = objfunc.Xmax(1);
    params.optimum =  0;

    global evaTime;
    evaTime = 0;
    SeedsMatrix         = zeros(params.seednum, params.dim);            % locations after individual movement
    fitness_best_array     = zeros(1,params.maxEva);                   % store the best fitness 

    % initialization
    for i = 1 : params.seednum
        SeedsMatrix(i,:) = repmat(params.lowerInit, 1, params.dim) + rand(1, params.dim).* repmat(params.upperInit - params.lowerInit, 1, params.dim);
    end

    % computing the fitness of each seeds
    SeedsFitness = zeros(1,params.seednum);

    SeedsFitness(1) = objfunc.obj_value(SeedsMatrix(1,:)');
    evaTime = evaTime + 1;
    fitness_best_array(evaTime) = SeedsFitness(1);

    for i=2:params.seednum
     SeedsFitness(i) = objfunc.obj_value(SeedsMatrix(i,:)');

     if(SeedsFitness(i)<fitness_best_array(evaTime))
        fitness_best_array(evaTime+1) = SeedsFitness(i);
     else
        fitness_best_array(evaTime+1) = fitness_best_array(evaTime);
     end
     evaTime = evaTime + 1;
    end



    % iteration by iteration
    while evaTime < params.maxEva

        % calculate the number of sons that each seed should generate
        sonsnum_array = sonsnum_cal(SeedsFitness,params);

        % calculate the exploding scope of sons  that each seed generate
        scope_array = scope_cal(SeedsFitness,params);

        % generate the sparks, based on the sparks number and explosion amplitude of each firework
        [SonsMatrix,SonsFitness,fitness_best_array] = sons_generate(sonsnum_array,scope_array,SeedsMatrix,params,fitness_best_array);

        % perform the gaussian mutation of seeds 
        [SeedsMatrixGauss,SeedsFitGaussMutation,fitness_best_array] = seedGaussMutation(SeedsMatrix, params,fitness_best_array);

        % all the seeds
        AllMatrix     = [SeedsMatrix;SonsMatrix;SeedsMatrixGauss]; 
        % attention SeedsFitness'
        AllFitness    = [SeedsFitness,SonsFitness,SeedsFitGaussMutation];

        % select the next iteration 
        [SeedsMatrix,SeedsFitnessCurrentIte]=selectNextIterationOnEntropy(AllMatrix,AllFitness,params);
        objfunc.draw_running(SeedsMatrix);

        %�ѽ��浽��һ����
        SeedsFitness = SeedsFitnessCurrentIte;
        if evaTime > params.maxEva
            evaTime = params.maxEva;
        end
    	[value,idx] = min(SeedsFitnessCurrentIte);
        outcome = [outcome; [evaTime value(1)]]; % �㷨���Ϊ���еľ���?
    end
    %% ****************==- �㷨�����ݴ��� -==*********************
    best_fvalue = value(1);
    best_solution = SeedsMatrix(idx(1));
    run_series = outcome;
    Altime = cputime - mytime;                  % ����ʱ��
    run_info = [Altime];
    
    fitness_best_array_return = fitness_best_array(1,params.maxEva_mod100);
    for i = 1 :params.maxEva_mod100
        fitness_best_array_return(i) = fitness_best_array(i*params.modStep);
    end
end

%% ����
function sonsnum_array = sonsnum_cal(fitness_array, params)

global Max_Sparks_Num;
global Min_Sparks_Num;
global Coef_Spark_Num;

fitness_max = max(fitness_array); 
fitness_sub_max = abs(fitness_max - fitness_array);
fitness_sub_max_sum = sum(fitness_sub_max);

sonsnum_array = zeros(1, params.seednum);

for i = 1 : params.seednum

    sonnum_temp = (fitness_sub_max(i) + eps) / (fitness_sub_max_sum + eps);
    sonnum_temp = round( sonnum_temp * Coef_Spark_Num);

    if sonnum_temp > Max_Sparks_Num
        sonnum_temp = Max_Sparks_Num;
    elseif sonnum_temp < Min_Sparks_Num
            sonnum_temp = Min_Sparks_Num;
    end
	
    sonsnum_array(i) = sonnum_temp; 
end
end

function scope_array = scope_cal(fitness_array, params)

global Coef_Explosion_Amplitude;

fitness_best = min(fitness_array);
fitness_sub_best = abs(fitness_best - fitness_array);
fitness_sub_best_sum = sum(fitness_sub_best);
scope_array = zeros(params.seednum);

for i=1:params.seednum
    scope_array(i) = Coef_Explosion_Amplitude * (fitness_sub_best(i) + eps) / (fitness_sub_best_sum + eps);  
end
end


function [sons_matrix, sons_fitness_array, fitness_best_array] = sons_generate(sonsnum_array, scope_array, seeds_matrix, params,fitness_best_array)

global total_fiteval_times;
global evaTime;

fiteval_time = sum(sonsnum_array);
total_fiteval_times = total_fiteval_times + fiteval_time;

sons_matrix = zeros(fiteval_time, params.dim);
sons_fitness_array = zeros(1, fiteval_time);    


sons_index = 1; 

for i = 1 : params.seednum
    for j = 1 : sonsnum_array(i)
	
        seed_position = seeds_matrix(i,:);
        
        allowed_dim_array = ones(1,params.dim);
        dimens_select = ceil(rand*params.dim); %select dimens_select
		
        offset = (rand*2-1) * scope_array(i); %Calculate the displacement:
        
        for k = 1 : dimens_select
            rand_dimen = ceil(rand*params.dim);
			
            while allowed_dim_array(rand_dimen)==0  
                rand_dimen = ceil(rand*params.dim);
            end
            allowed_dim_array(rand_dimen)=0;
			
            seed_position(rand_dimen) = seed_position(rand_dimen) + offset;
            
            if seed_position(rand_dimen) > params.upperBound || seed_position(rand_dimen) < params.lowerBound
			% map to the search range of the sparks which are out of the bound
                span = params.upperBound - params.lowerBound;
                seed_position(rand_dimen) = params.lowerBound + rem(abs(seed_position(rand_dimen)), span);
            end
        end

        
        sons_matrix(sons_index,:) = seed_position;
        sons_index = sons_index + 1;
    end
end

for i = 1 : fiteval_time
    global objfunc ;
    sons_fitness_array(i) = objfunc.obj_value(sons_matrix(i,:)');
	 if(sons_fitness_array(i)<fitness_best_array(evaTime))
		fitness_best_array(evaTime+1) = sons_fitness_array(i);
	 else
		fitness_best_array(evaTime+1) = fitness_best_array(evaTime);
	 end
	 evaTime = evaTime + 1;
end
end

function [seed_gaussian_matrix,fit_seed_gaussian_array,fitness_best_array] = seedGaussMutation(seedMatrix, params,fitness_best_array)

    gaussian_number = params.gaussianNum;
    global total_fiteval_times;
    global objfunc ;
    global evaTime;

    total_fiteval_times = total_fiteval_times + gaussian_number;

    dim = params.dim;

    seed_gaussian_matrix = zeros(gaussian_number, dim);
    fit_seed_gaussian_array=zeros(1, gaussian_number);


     for i = 1 : gaussian_number
            % select a firework
            seed_select = ceil(rand * params.seednum);
            rand_dimens = ceil(rand * dim);

            seed_gaussian_matrix(i,:)= seedMatrix(seed_select,:);
            coef_gaussian         = normrnd(0,1)+1;
            allow_dimension_array     = ones(1,dim);

            for j = 1:rand_dimens

                dim_mutation = ceil(rand * dim);
                while allow_dimension_array(dim_mutation)==0
                    dim_mutation = ceil(rand * dim);
                end
                allow_dimension_array(dim_mutation) = 0;
                seed_gaussian_matrix(i, dim_mutation) = seed_gaussian_matrix(i, dim_mutation) * coef_gaussian;

                % map to the search range of the sparks which are out of the bound
                if seed_gaussian_matrix(i, dim_mutation) > params.upperBound || seed_gaussian_matrix(i, dim_mutation) < params.lowerBound
                    span = params.upperBound - params.lowerBound;
                    seed_gaussian_matrix(i,dim_mutation) = params.lowerBound + rem(abs(seed_gaussian_matrix(i, dim_mutation)), span);
                end
            end

            fit_seed_gaussian_array(i) =  objfunc.obj_value(seed_gaussian_matrix(i,:)');  
            if(fit_seed_gaussian_array(i)<fitness_best_array(evaTime))
                fitness_best_array(evaTime+1) = fit_seed_gaussian_array(i);
             else
                fitness_best_array(evaTime+1) = fitness_best_array(evaTime);
             end
             evaTime = evaTime + 1;
     end
end

 function [seedNextMatrix,seedNextMatrixFitness]=selectNextIterationOnEntropy(seedMatrix,seedFitness,params)
    %  		select the fireworks for next iteration
    %       


    particleNum = size(seedFitness,2);
    dim = size(seedMatrix,2);

    [~,tempIndex]=sort(seedFitness);
    bestIndex = tempIndex(1);

    %  Prelocate the parameters

    seedNextMatrix              = zeros(params.seednum,dim);%�������?
    seedNextMatrixFitness       = zeros(1,params.seednum);    

    NormMatrix                  = zeros(particleNum,particleNum);
    DistanceArry                = zeros(1,particleNum);
    ProbabiSelctPaticle         = zeros(1,particleNum);
    addProbabiSelctPaticle      = zeros(1,particleNum);

    randPoint                   = rand(1,params.seednum-1);
    randPointIndex              = zeros(1,params.seednum-1);


    for i=1:particleNum
       for j=i:particleNum
           % compute the distace between each two
           NormMatrix(i,j)=mathNorm(seedMatrix(i,:),seedMatrix(j,:));
       end
    end
    NormMatrix = NormMatrix + NormMatrix';

    % compute the sum distance for a spark to others 
    for i=1:particleNum
        for j=1:particleNum
            DistanceArry(i)=DistanceArry(i)+NormMatrix(i,j);
        end
    end

    % probobility computation
    for i=1:particleNum
        ProbabiSelctPaticle(i)=DistanceArry(i)/sum(DistanceArry);
    end



    %selection perform
    addProbabiSelctPaticle(1)=ProbabiSelctPaticle(1);
    for i=2:particleNum
        addProbabiSelctPaticle(i)=ProbabiSelctPaticle(i)+addProbabiSelctPaticle(i-1);
    end


    for i=1:params.seednum-1
        j=1;
        while(randPoint(i)>addProbabiSelctPaticle(j) && j<particleNum)
            j=j+1;
        end
        randPointIndex(i)=j;
    end

    % the best
    seedNextMatrix(1,:)=seedMatrix(bestIndex,:);
    seedNextMatrixFitness(1)=seedFitness(bestIndex);

    % the others
    for i = 1:params.seednum - 1
        seedNextMatrix(i+1,:) = seedMatrix(randPointIndex(i),:);
        seedNextMatrixFitness(i+1) = seedFitness(randPointIndex(i));
    end
 end

function NormDistance=mathNorm(arrayX,arrayY)
    % ˵��������������ֻ�ܼ����������ķ�ʽ���룬û�������к��е��ж�

    %%
    xsize=size(arrayX,2);
    ysize=size(arrayY,2);
    NormDistance2=0;


    %%


    if xsize==ysize
        for i=1:xsize
           NormDistance2=NormDistance2+ (arrayX(i)-arrayY(i))*(arrayX(i)-arrayY(i));
        end
        NormDistance=sqrt(NormDistance2);
    else
       error('Computing mathNorm of arrayX and arrayY : the dimention is not match. ') ;
    end
end