function [best_fvalue,best_solution,run_series,run_info,user_set] = SE(func_name , algo_config,func_config,user_config)

%**************************************************************************
% 
% 
% �˴�дע��
% 
%**************************************************************************
    format long;
    format compact;
    rng('shuffle'); 
    mytime = cputime;            % ��ʼʱ��
    %% *****************==- ��ʼ�����Ժ��� -==**********************
    if exist('func_config','var')
        objfunc = feval(func_name,func_config);  % ʵ��һ��Ŀ�꺯����
    else
        objfunc = feval(func_name);  % ʵ��һ��Ŀ�꺯����
        func_config = [];
    end
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
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
    globalBestValue = Inf;
    localBestValue = Inf;
    swarm = objfunc.initialize( PopuSize );  % ��ʼ����Ⱥ,���һ����Ⱥ��ģΪPopuSize����Ⱥ 
    f_swarm = objfunc.obj_value( swarm );
    newswarm = swarm;

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
        
        if (fitval(1,1) <= globalBestValue)
            globalBestValue = fitval(1,1);            % ��¼ȫ������ֵ?
            globalBestSolution = swarm(:,fitindex(1));  % ��¼ȫ�����Ž�
        end
        
        FES = FES + PopuSize;                       % �����������?
        outcome = [outcome; [FES globalBestValue]]; % �㷨���Ϊ���еľ���?
        if FES >= maxFES
          break;          % �����������������?
        end
%         objfunc.draw_running(swarm,'PSO');

        
        
       %% ��ɢ����
            for i = 1:PopuSize
                r1 = rand();
                r2 = 1 + r1;
                r3 = randn(DIM,1);
                r4 = 1 - r3;
                r5 = rand(DIM,1);
%                     newswarm(:,i) = -r5.*swarm(:,i) + (1+r5).*((r3.*swarm(:,randi(PopuSize)) + r4.*globalBestSolution(:,1)));
%                     newswarm(:,i) = ((1-2*r5).*swarm(:,i) + 2*r5.*((r3.*localBestSolution(:,1) + r4.*globalBestSolution(:,1))));
                    newswarm(:,i) = ((1-2*r5).*swarm(:,i) + 2*r5.*((r3.*swarm(:,randi(PopuSize)) + r4.*globalBestSolution(:,1))));
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

