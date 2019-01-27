function [best_fvalue,best_solution,run_series,run_info,user_set] = IGHS(func_name , algo_config,func_config,user_config)
%****************************************************************************
%Reference:  
    % A. K. Qin, V. L. Huang, and P. N. Suganthan,��Differential evolution
    % algorithm with strategy adaptation for global numerical optimization,��
    % IEEE Trans. Evolut. Comput., vol. 13, no. 2, pp. 398�C417, Apr. 2009.
% Note: 
    % We obtained the MATLAB source code from the authors, and did some
    % minor revisions in order to solve the 25 benchmark test functions,
    % however, the main body was not changed.
%*****************************************************************************
    format long;
    format compact;
    rand('seed', sum(100*clock));% ���������
    getrand = @(l, r) l + (r - l) .* rand();%%�ڼ������ȡ
    getnum = @(num) fix(1 + rand() * num);%%�ڼ������ȡ
    mytime = cputime;
    outcome = [];

    %% *****************==- ��ʼ�����Ժ��� -==**********************
    if exist('func_config','var')
        objfunc = feval(func_name,func_config);  % ʵ��һ��Ŀ�꺯����
    else
        objfunc = feval(func_name);  % ʵ��һ��Ŀ�꺯����
        func_config = [];
    end
    objfunc = feval(func_name,func_config);
    N = objfunc.Xdim;   %%�����ά��
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
    
    %% *****************==- ��ʼ���㷨���� -==**********************
    if ~exist('algo_config','var')
        algo_config = [];  % �����ʼ��
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

    HMCR = 0.995;%%�������ѡȡ�ĸ���
    PAR = 0.4;%%΢���ĸ���

    ansr = zeros(HMS, N);%%����⣬HMS*N��С�ľ���
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
	range = [Xmin, Xmax];%%��ֵ��ȡֵ��Χ
    F = zeros(1, HMS);%%��¼Ŀ�꺯���С

    %%��ʼ�������
    for i = 1 : HMS
        for j = 1 : N
            ansr(i, j) = getrand(range(1), range(2));%%����ʼ��ÿ������
        end
        F(i) = [objfunc.obj_value( ansr(i, :)')]';%%��¼ÿ��������Ŀ�꺯���ֵ
    end
    
    [fitval, fitindex] = min(F);
    globalBestValue = fitval(1,1);   
    globalBestSolution = ansr(fitindex(1),:);
    outcome = [HMS globalBestValue];

    ANSP = zeros(NI, N);    %%��¼������в����ÿ������
    FP = zeros(1, NI);      %%��¼������в����ÿ��������Ŀ�꺯��ֵ
    F_best = zeros(1 , NI); %%�洢ÿ�ε���е����Ž�

    %% *****************==- ��ѭ�� -==**********************
    for T = 1 : NI
        new_solution = zeros(1, N); %% ��ʱ��¼���������

        %�ҳ�Ŀǰ������е���������Ž�
        worst_idx = 1;  %%worst_idx ��¼��������
        best_idx = 1;   %%best_idx  ��¼���Ž�����
        for i = 1 : HMS
            if (F(worst_idx) < F(i))
                worst_idx = i;
            end
            if (F(best_idx) > F(i))
                best_idx = i;
            end
        end

        %����µĽ�����
        for i = 1 : N
            pos = rand();
            if (pos < HMCR)
                 new_solution(i) = 2 * ansr(best_idx , i) - ansr(worst_idx , i);
                 if(new_solution(i) > range(2))
                     new_solution(i) = range(2);
                 else if(new_solution(i) < range(1))
                         new_solution(i) = range(1);
                     end
                 end
                 new_solution(i) = ansr(worst_idx , i) + rand() * (new_solution(i) - ansr(worst_idx , i));%%�ڿ������һ���µ�����
                 posr = rand();
                 if (posr < PAR)
                     new_solution(i) = ansr(best_idx , getnum(N));%%�����ŵ�һ�������������ѡ��
                 end
            else
                new_solution(i) = getrand(range(1), range(2));%%�ڿ���������
            end
        end

        %�жϸ��¼����
        F_new = [objfunc.obj_value( new_solution')]'; %%�µ�������Ŀ�꺯��ֵ
        if (F_new < F(worst_idx))
            F(worst_idx) = F_new;  %%Ŀ�꺯��ֵ����
            ansr(worst_idx, :) = new_solution;%%��������
        end
        
		if (globalBestValue >= F_new)
            globalBestValue = F_new;
            globalBestSolution(1,:) = new_solution;
        end
        if T > NI
            break;
        end        
        if(mod(T,100) == 0)
            outcome = [outcome ; (T) globalBestValue];
        end ;
        

       %�洢�������ɵ����� 
       ANSP(T, :) = new_solution;
       FP(T) = F_new;
       F_best(T) = F(best_idx);%%�洢ÿ�ε���е����Ž�
    end
    
    %% ****************==- �㷨�����ݴ��� -==*********************
    best_fvalue = globalBestValue;
    best_solution = globalBestSolution;
    run_series = outcome;
    Altime = cputime - mytime;                  % ����ʱ��
    run_info = [Altime];
end
    