function outcome = DSRSL(func_name,PopuSize,maxFES)
%**************************************************************************
% 
% 
% 此处写注释
% 
%**************************************************************************
    format long;
    format compact;
    rand('seed', sum(100*clock));% 随机数种子
    mytime = cputime;            % 初始时间
    switch nargin
        case 1
            PopuSize = 40;
            maxFES = 60000;
        case 2
            maxFES = 60000;
    end
    FES = 0;
    outcome = [];
    objfunc = feval(func_name);  % 实例化一个目标函数类,这样就可以调用这个类的函数了
    DIM = objfunc.Xdim;          % 获取目标函数类的一些参数
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
    globalBestValue = Inf;
    localBestValue = Inf;
    swarm = objfunc.argument( PopuSize );  % 初始化种群,生成一个种群规模为PopuSize的种群 
    newswarm = swarm;
    document = swarm;
    doc_value = inf*ones(1,PopuSize);
    f_swarm = inf*ones(1,PopuSize);
    change = zeros(1,PopuSize);
    unchange = 60; 
    gagen = 0;
    simgae = 0;
    loc_change = zeros(1,PopuSize);
%     loc_change = 0;
    loc_change_num = 0;
    doc_change = 0;
    doc_change_num = 0;
    % 主循环
    while FES < maxFES
        f_swarm_new = objfunc.fun_fitness( newswarm ); % 计算目标函数值,注意这里是函数值并非适应度值

       %% 选择
        for i = 1:PopuSize
            if(f_swarm(i) >= f_swarm_new(i))
               swarm(:,i) = newswarm(:,i);    
               f_swarm(i) = f_swarm_new(i);
            end
        end
        [fitval, fitindex] = min(f_swarm);
        
        gen = mod(FES/PopuSize,PopuSize) + 1;
        loc_change(gen) = 0;      
        if (fitval(1,1) <= globalBestValue)
            globalBestValue = fitval(1,1);            % 记录全局最优值?
            globalBestSolution = swarm(:,fitindex(1));  % 记录全局最优解
        end
        if (fitval(1,1) <= localBestValue)
            localBestValue = fitval(1,1);            % 记录当前种群最优值?
            localBestSolution = swarm(:,fitindex(1));  % 记录当前种群最优解
            loc_change(gen) = loc_change(gen) + 1;
            loc_change_num = loc_change_num + 1;
        end
        
        FES = FES + PopuSize;                       % 更新运算次数
        outcome = [outcome; [FES globalBestValue]]; % 算法输出为两列的矩阵
        if FES >= maxFES
          break;          % 若超出最大迭代次数跳出
        end
%         %% 绘图
%         hold off;
%         scatter(swarm(1,:),swarm(2,:));
%         hold on;
%         objfunc.fun_figue(3);   
%         pause(.001);
%         
        %% 存档
        gen = mod(FES/PopuSize,PopuSize) + 1;
        change(gen) = 0;
        doc_change_num = doc_change_num + 1;
        doc_change = mod(doc_change_num,PopuSize) + 1;
        if( doc_value(doc_change) > globalBestValue)
           document(:,doc_change) = globalBestSolution(:,1);    
           doc_value(doc_change) = globalBestValue;
           change(gen) = change(gen) + 1;
           doc_change = doc_change + 1;
        else
            doc_change_num = doc_change_num - 1;
        end
%         gen = mod(FES/PopuSize,PopuSize) + 1;
%         change(gen) = 0;
%         doc_change_num = doc_change_num + 1;
%         doc_change = mod(doc_change_num,PopuSize) + 1;
%         if( doc_value(doc_change) > localBestValue)
%            document(:,doc_change) = localBestSolution(:,1);    
%            doc_value(doc_change) = localBestValue;
%            change(gen) = change(gen) + 1;
%            doc_change = doc_change + 1;
%         else
%             doc_change_num = doc_change_num - 1;
%         end
        
        %% 重置新的种群 && (PopuSize/10sum(loc_change(mod(loc_change_num,PopuSize) + 1)) < 1)
        if ((sum(loc_change) < 1)  && (FES > PopuSize)) || rand < 0.3
            if(doc_change_num < PopuSize)
%                 swarm( :,1:floor(PopuSize/2) ) = objfunc.argument( floor(PopuSize/2) );
                swarm( :,1:PopuSize) = objfunc.argument( PopuSize );
                f_swarm = objfunc.fun_fitness( swarm );
            end
            change = ones(1,PopuSize);
            localBestValue = Inf;
        else
%             for i = 1:PopuSize
%                swarm(:,i) = document(:,i);    
%                f_swarm(i) = doc_value(i);  
%                unchange = unchange - 1;
%             end
               swarm = document;    
               f_swarm = doc_value;  
            if rand < 0.6
       
            else
                swarm( : , 1:floor(PopuSize/2) ) = objfunc.argument( floor(PopuSize/2) );
                f_swarm = objfunc.fun_fitness( swarm ); 
            end
            change = ones(1,PopuSize);
            localBestValue = Inf;
        end
            simgae = simgae + 1;
            %% 单纯形代
            for i = 1:PopuSize
                if i < PopuSize - 1
                    index(1) = i;
                    index(2) = i + 1;
                    index(3) = i + 2;
                else
                    index(1) = i;
                    index(2) = mod(mod(i + 1,PopuSize),PopuSize + 1) + 1;
                    index(3) = mod(mod(i + 2,PopuSize),PopuSize + 1) + 1;

                end
                [val_max,index_max] = max([f_swarm(index(1)) f_swarm(index(2)) f_swarm(index(3))]);
                [val_min,index_min] = min([f_swarm(index(1)) f_swarm(index(2)) f_swarm(index(3))]);                
                equal = 0;
                if index_max == index_min
                    a = index(1);
                    c = index(2);
                    b = index(3); 
                    equal = 1;
                else
                    a = index(index_max);
                    c = index(index_min);
                    index([index_min index_max]) = [];
                    b = index;                   
                end
                
                r1 = rand();
                r2 = rand()*0.3;
                r3 = 1 - r2;
                if rand < 0.8 &&  equal == 1
                    newswarm(:,i) = swarm(:,c) + r1*2*((r2*swarm(:,b) + r3*swarm(:,a))/2 - swarm(:,c));
%                 else
                elseif rand < 0.6
                      newswarm(:,i) = (swarm(:,i) + r1*2*((localBestSolution(:,1) + globalBestSolution(:,1))/2 - swarm(:,i)));
%                     newswarm(:,i) = posCross'.*document(:,i) + posCross_'.*(swarm(:,i) + r1*2*((localBestSolution(:,1) + globalBestSolution(:,1))/2 - swarm(:,i)));
                else
                    newswarm(:,i) = swarm(:,a) + r1*2*((swarm(:,b) + swarm(:,c))/2 - swarm(:,a));
                end
                
                rd = randi( DIM,1,2 );
                newswarm(rd,i) = globalBestSolution(rd,1);
%                 for mi = 1:PopuSize
%                     posCross = (rand(1,DIM)<0.3);
%                     posCross_ = 1 - posCross;
%                     newswarm(:,mi) = posCross'.*document(:,mi) + posCross_'.*newswarm(:,mi);
%                 end
            end
            indexLB = find( newswarm < Xmin );
            newswarm(indexLB) = min( Xmax, 2*Xmin - newswarm(indexLB) );
            indexUB = find( newswarm > Xmax );
            newswarm(indexUB) = max( Xmin, 2*Xmax - newswarm(indexUB) );
    end
save localBestSolution localBestSolution 
    Altime = cputime - mytime;                  % 计算时间
    display(['      DSRSL time = ' num2str(Altime)]);  % 显示时间
    disp(['      最优值  =  ' num2str(globalBestValue,'%10.5e')  ]);
end