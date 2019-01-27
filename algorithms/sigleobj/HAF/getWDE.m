function WDE = getWDE(swarm,PopuSize,DIM,globalBestSolution)
        F = normrnd(0.5, 0.3, 1, PopuSize); %0.5 + rand*0.1;
        rr = (rand(1, DIM)>0.5);  %different CR with 0.9 or 0.2
        CR = 0.9.*rr + 0.2.*(1-rr); %0.9 + rand * 0.05;
        for mi = 1:PopuSize
            rIn = RandIndex(PopuSize,2);
            Vswarm(:,mi) = globalBestSolution+F(mi)*(swarm(:,rIn(1))-swarm(:,rIn(2)));
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
        WDE = swarm\Uswarm;
end

function randindex = RandIndex(pSize,k)
%size: input popusize;  k: k different subindex are required
%randindex: k different subindex are returned in this array
randpert = randperm(pSize);
randindex = randpert(1:k);
end
