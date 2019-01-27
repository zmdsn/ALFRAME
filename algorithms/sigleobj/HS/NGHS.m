function [best_fvalue,best_solution,run_series,run_info,user_set] = NGHS(func_name , algo_config,func_config,user_config)
    getrand = @(l, r) l + (r - l) .* rand();%%�ڼ�������?
    mytime = cputime;
    %% *****************==- ��ʼ�����Ժ��� -==**********************
    if exist('func_config','var')
        objfunc = feval(func_name,func_config);  % ʵ��һ��Ŀ�꺯����
    else
        objfunc = feval(func_name);  % ʵ��һ��Ŀ�꺯����
        func_config = [];
    end
    objfunc = feval(func_name,func_config);
    N = objfunc.Xdim;   %%�����ά��?
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
    
    %% *****************==- ��ʼ���㷨���� -==**********************
    if ~exist('algo_config','var')
        algo_config = [];  % �����ʼ��?
    end
    % ��Ⱥ��ģ
    if isfield(algo_config,'PopuSize')
        HMS = algo_config.PopuSize;
    else
        HMS = 5;
    end
    % ������
    if isfield(algo_config,'maxFES')
        NI = algo_config.maxFES;
    else
        NI = 60000;
    end
    Pm = 0.005; %%����ĸ���?

    ansr = zeros(HMS, N);%%����⣬HMS*N��С�ľ���
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
	range = [Xmin, Xmax];%%��ֵ��ȡֵ��Χ
    F = zeros(1, HMS);%%��¼Ŀ�꺯����?

    %%��ʼ�������?
    for i = 1 : HMS
        for j = 1 : N
            ansr(i, j) = getrand(range(1), range(2));%%����ʼ��ÿ������
        end
        F(i) = objfunc.obj_value( ansr(i, :)');%%��¼ÿ��������Ŀ�꺯����?
    end
    [fitval, fitindex] = min(F);
    globalBestValue = fitval(1,1);   
    globalBestSolution = ansr(fitindex(1),:);
    outcome = [HMS globalBestValue];

%     ANSP = zeros(NI, N);%%��¼������в����ÿ������
%     FP = zeros(1, NI);%%��¼������в����ÿ��������Ŀ�꺯��ֵ

    %% *****************==- Main Loop -==**********************
    for T = 1 : NI
        ans_temp = zeros(1, N);%%��ʱ��¼���������?
        %%�ҳ�Ŀǰ������е���������Ž�
        worst = 1;%%num_i1��¼��������
        best = 1;%%num_i2��¼���Ž�����
        for i = 1 : HMS
            if (F(i) > F(worst))
                worst = i;
            elseif (F(i) < F(best))
                best = i;
            end
        end
%         [v ,worst] = max(F);
%         [v ,best] = min(F);
        %%����µ�����?
        for i = 1 : N
            ans_temp(i) = 2* ansr(best , i) - ansr(worst , i);
            if(ans_temp(i) > range(2))
                ans_temp(i) = range(2);
            elseif(ans_temp(i) < range(1))
                ans_temp(i) = range(1);
            end
            r = rand*2;
            
            ans_temp(i) = ansr(worst , i) + rand * ( ans_temp(i) - ansr(worst , i) );%%�ڿ������һ���µ�����?
%             ans_temp(i) =  r * ansr(best , i) + (1-r) * ansr(worst , i) ;%%�ڿ������һ���µ�����?
%             ans_temp(i) =  rand *2* ansr(best , i) + (1-rand*2) * ansr(worst , i) ;%%�ڿ������һ���µ�����?
            posr = rand();
            if (posr < Pm)
                ans_temp(i) = Xmin + rand*(Xmax - Xmin);%%�ڿ���������
            end
        end

        %%�жϸ��¼����?
        F_temp = objfunc.obj_value( ans_temp');%%�µ�������Ŀ�꺯��ֵ
%         if (F_temp < F(worst))
            F(worst) = F_temp;%%Ŀ�꺯��ֵ����
            ansr(worst, :) = ans_temp;%%��������
%         end
        
        if globalBestValue > F_temp
            globalBestValue = F_temp;
            globalBestSolution(1,:) = ans_temp;
        end
        if T + HMS >= NI
            T = NI - HMS ;
            break;          % �����������������?
        end
        globalBestValue
        if(mod(T,100) == 0)
            outcome = [outcome ; (T + HMS) globalBestValue];
        end 
%         objfunc.draw_running( ansr' );
   end
    

    % plot3(ANSP(1 : 10 : NI, 1), ANSP(1 : 10 : NI, 2), ANSP(1 : 10 : NI, 3), 'r.', 'markersize', 1);

    %% ****************==- �㷨�����ݴ��� -==*********************
    best_fvalue = globalBestValue;
    best_solution = globalBestSolution;
    run_series = outcome;
    Altime = cputime - mytime;                  % ����ʱ��
    run_info = [Altime];
end
