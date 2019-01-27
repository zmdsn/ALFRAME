function [best_fvalue,best_solution,run_series,run_info,user_set] = GAD(func_name , algo_config,func_config,user_config)

%**************************************************************************
% 
% 
% �˴�дע��
% 
%**************************************************************************
    format long;
    format compact;
    rand('seed', sum(100*clock));% ���������?
    mytime = cputime;            % ��ʼʱ��
    %% *****************==- ��ʼ�����Ժ��� -==**********************
    if exist('func_config','var')
        objfunc = feval(func_name,func_config);  % ʵ��һ��Ŀ�꺯����
    else
        objfunc = feval(func_name);  % ʵ��һ��Ŀ�꺯����
        func_config = [];
    end
    objfunc = feval(func_name,func_config);
    DIM = objfunc.Xdim;          % ��ȡĿ�꺯�����һЩ����?
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
    
    %% *****************==- ��ʼ���㷨���� -==**********************
    if ~exist('algo_config','var')
        algo_config = [];  % �����ʼ��?
    end
    % ��Ⱥ��ģ
    if isfield(algo_config,'PopuSize')
        PopuSize = algo_config.PopuSize;
    else
        PopuSize = 50;
    end
    % ������
    if isfield(algo_config,'maxFES')
        maxFES = algo_config.maxFES;
    else
        maxFES = 600000;
    end
    
    FES = PopuSize;
    outcome = [];
    objfunc = feval(func_name);  % ʵ��һ��Ŀ�꺯����,����Ϳ��Ե��������ĺ�����?
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
    globalBestValue = Inf;
    localBestValue = Inf;
    swarm = objfunc.initialize( PopuSize );  % ��ʼ����Ⱥ,���һ����Ⱥ��ģΪPopuSize����Ⱥ 
    f_swarm = objfunc.obj_value( swarm );
    newswarm = swarm;
    document = swarm;
    doc_value = inf*ones(1,PopuSize);
    doc_change = 1;
    change = zeros(1,PopuSize);
    while FES < maxFES
        f_swarm_new = objfunc.obj_value( newswarm ); % ����Ŀ�꺯��ֵ,ע�������Ǻ���ֵ������Ӧ��ֵ
       %% ѡ��
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
            globalBestValue = fitval(1,1);            % ��¼ȫ������ֵ?
            globalBestSolution = swarm(:,fitindex(1));  % ��¼ȫ�����Ž�
        end
        if (fitval(1,1) <= localBestValue)
            localBestValue = fitval(1,1);            % ��¼��ǰ��Ⱥ����ֵ?
            localBestSolution = swarm(:,fitindex(1));  % ��¼��ǰ��Ⱥ���Ž�
            loc_change(gen) = loc_change(gen) + 1;
        end
        FES = FES + PopuSize;                       % �����������?
        outcome = [outcome; [FES globalBestValue]]; % �㷨���Ϊ���еľ���?
        if FES >= maxFES
          break;          % �����������������?
        end
% %         %% ��ͼ
%         if nargin >= 4 && draw == 1
%             hold off;
%             scatter(swarm(1,:),swarm(2,:));
%             hold on;
%             scatter(document(1,:),document(2,:));
%             hold on;
%             objfunc.fun_figue(3);   
%             pause(.001);
%         end
%         
        %% �浵
        gen = mod(FES/PopuSize,PopuSize) + 1;
        change(gen) = 0;
        if( doc_value(doc_change) > globalBestValue)
           document(:,doc_change) = globalBestSolution(:,1);    
           doc_value(doc_change) = globalBestValue;
           change(gen) = change(gen) + 1;
           doc_change = doc_change + 1;
           if doc_change > PopuSize
               doc_change = 1;
           end
        end     
        %% �����µ���Ⱥ
        if ((sum(loc_change) <= PopuSize/10))
            swarm = document;    
            f_swarm = doc_value;                
            change = ones(1,PopuSize);     
            localBestValue = Inf;
        elseif  rand < 1/PopuSize % rand < 0.1
            if rand < 0.5
                swarm = document;
                f_swarm = doc_value;
                change = ones(1,PopuSize);
                localBestValue = Inf;
            elseif rand < 0.5
                swarm = document + rand*PopuSize/FES*objfunc.initialize( PopuSize );
                f_swarm = objfunc.obj_value( swarm );  
                FES = FES + PopuSize;  
                outcome = [outcome; [FES globalBestValue]]; 
                change = ones(1,PopuSize);
                localBestValue = Inf;
            else 
                swarm = objfunc.initialize( PopuSize );
                f_swarm = objfunc.obj_value( swarm ); 
                FES = FES + PopuSize;  
                outcome = [outcome; [FES globalBestValue]]; 
                change = ones(1,PopuSize);
                localBestValue = Inf;
            end
        end
       %% ��ɢ����
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
                    worst = index(1);
                    best = index(2);
                    middle = index(3); 
                    equal = 1;
                else
                    worst = index(index_max);
                    best = index(index_min);
                    index([index_min index_max]) = [];
                    middle = index;                   
                end
                
%                 r1 = rand();
%                 r2 = rand();
%                 r3 = 1 - r2;
% %                 if equal == 1
%                     newswarm(:,i) = swarm(:,best) + r1*2*((r2*swarm(:,middle) + r3*swarm(:,worst))/2 - swarm(:,best));
%                 elseif rand < 0.5
%                     newswarm(:,i) = (swarm(:,i) + r1*2*((localBestSolution(:,1) + globalBestSolution(:,1))/2 - swarm(:,i)));
%                 else
%                     newswarm(:,i) = swarm(:,worst) + r1*2*((swarm(:,middle) + swarm(:,best))/2 - swarm(:,worst));
%                 end

                r1 = rand();
                r2 = 1 + r1;
                r3 = rand();
                r4 = 1 - r3;
%                 r5 = rand()*0.5+0.5;
                r5 = rand();
%                 if equal == 1
%                     newswarm(:,i) = r1*swarm(:,i) + (1-2*r1)*((globalBestSolution(:,1)));
% %                     newswarm(:,i) = r1*swarm(:,best) + (1-2*r1)*((swarm(:,middle) + swarm(:,worst)));
% %                     newswarm(:,i) = r1*swarm(:,best) + (1-2*r1)*((r3*swarm(:,middle) + r4*swarm(:,worst)));
%                 elseif rand < 0.5
%                     newswarm(:,i) = (r5*swarm(:,i) + (1-2*r5)*((localBestSolution(:,1) + globalBestSolution(:,1))));
%                     newswarm(:,i) = (swarm(:,i) + r1*2*((localBestSolution(:,1) + globalBestSolution(:,1))/2 - swarm(:,i)));
%                     newswarm(:,i) = r1*swarm(:,best) + 2*r2*PopuSize/FES*((r3*localBestSolution(:,1) + r4*globalBestSolution(:,1)));
%                 else
% %                     newswarm(:,i) = (1-2*r5)*swarm(:,worst) + r5*((r3*swarm(:,middle) + r4*swarm(:,best)));
% %                     newswarm(:,i) = (1-2*r5)*swarm(:,worst) + r5*((swarm(:,middle) + swarm(:,best)));
%                     newswarm(:,i) = ((1-2*r5)*swarm(:,i) + r5*((localBestSolution(:,1) + globalBestSolution(:,1))));
%                 end
%                 if equal == 1
%                     newswarm(:,i) = r1*swarm(:,i) + (1-2*r1)*((globalBestSolution(:,1)));
%                 else
%                     newswarm(:,i) = ((1-2*r5)*swarm(:,i) + r5*((localBestSolution(:,1) + globalBestSolution(:,1))));
%                 end
                if    equal == 1   
%                     newswarm(:,i) = 2*r1*swarm(:,i) + (1-2*r1)*((globalBestSolution(:,1)));
%                     newswarm(:,i) = 2*r1*swarm(:,i) + (1-2*r1)*((globalBestSolution(:,1)));
                    newswarm(:,i) = ((1-2*r5)*swarm(:,best) +  r5*((r3*swarm(:,worst) + r4*swarm(:,middle))));
                else
                    newswarm(:,i) = ((1-2*r5).*swarm(:,i) + 2*r5.*((r3*localBestSolution(:,1) + r4*globalBestSolution(:,1))));

%                     newswarm(:,i) = ((1-2*r5)*swarm(:,i) + 2*r5*((r3*localBestSolution(:,1) + r4*globalBestSolution(:,1))));
                end
            end
            % Խ�紦��
            indexLB = find( newswarm < Xmin );
            newswarm(indexLB) = min( Xmax, 2*Xmin - newswarm(indexLB) );
            indexUB = find( newswarm > Xmax );
            newswarm(indexUB) = max( Xmin, 2*Xmax - newswarm(indexUB) );
    end
   %% ****************==- �㷨�����ݴ��� -==*********************
    best_fvalue = globalBestValue;
    best_solution = globalBestSolution;
    run_series = outcome;
    Altime = cputime - mytime;                  % ����ʱ��
    run_info = [Altime];    
    
%     Altime = cputime - mytime;                  % ����ʱ��
%     display(['      GAD time = ' num2str(Altime)]);  % ��ʾʱ��
%     disp(['      ����ֵ  =  ' num2str(globalBestValue,'%10.5e')  ]);
end

