function WGA = getWGA(swarm,f_swarm,PopuSize,DIM)
    pCrossover = 0.8;  
    pMuation = 0.1; 
    % Selection
    % roulette gambling
    ifs = 1./( f_swarm ) ;
    sumfitness = sum( ifs );
    sumf = ifs ./sumfitness;
    index=[];
    for i=1:PopuSize   %转sizepop次轮盘
        pick=rand;
        while pick==0
            pick=rand;
        end
        for j=1:PopuSize
            pick=pick-sumf(j);
            if pick< 0
                index=[index j];
                break;  %寻找落入的区间，此次转轮盘选中了染色体i，注意：在转sizepop次轮盘的过程中，有可能会重复选择某些染色体
            end
        end
    end
    parent = swarm(:,index);

    % Crossover
    for mi = 1:PopuSize
        randomIndex = randperm(PopuSize);
        % select two solutions for crossover:
        if rand < pCrossover
           ii1 = randomIndex(1); 
           ii2 = randomIndex(2);
           weight = rand(DIM,1)< 0.5; 
           newswarm(:,mi) = weight.*parent(:,ii1) + (1-weight).*parent(:,ii2); 
        else
           newswarm(:,mi) = swarm(:,mi);
        end
    end

    % Mutation:
    for mi = 1 : PopuSize
        for nj =1: DIM
            if rand < pMuation
               newswarm(nj,mi) = newswarm(nj,mi) + randn;
            end
        end
    end
    
    WGA = swarm\newswarm;
end