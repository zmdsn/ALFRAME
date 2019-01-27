function [WPSO,v] = getWPSO(swarm,v,bestL,globalBestSolution,PopuSize,DIM)
    w = 1/(2*log(2));
    c1 = 0.5 + log(2);
    c2 = c1;
    for i=1:PopuSize
        v(:,i) = w*v(:,i) + c1*rand(DIM,1).*(bestL(:,i)-swarm(:,i)) + c2*rand(DIM,1).*(globalBestSolution(:,1)-swarm(:,i));
        newPAR(:,i) = swarm(:,i) + v(:,i);
    end      
    WPSO = swarm\newPAR;
end
