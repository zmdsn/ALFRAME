function [best_fvalue,best_solution,run_series,run_info,user_set] = IGHS(func_name , algo_config,func_config,user_config)
    mytime = cputime;
    getrand = @(l, r) l + (r - l) .* rand();%%�ڼ�������?
    getnum = @(num) fix(1 + rand() * num);%%�ڼ�������?
    dt = @(num) 2 * num * (rand() - 0.5);%%�������΢��?

    %% *****************==- ��ʼ�����Ժ��� -==**********************
    if exist('func_config','var')
        objfunc = feval(func_name,func_config);  % ʵ��һ��Ŀ�꺯����
    else
        objfunc = feval(func_name);  % ʵ��һ��Ŀ�꺯����
        func_config = [];
    end
    objfunc = feval(func_name,func_config);
    DIM = objfunc.Xdim;   %%�����ά��?
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
        NI = 40;
    end
    outcome = [];
    
    HMCR = 0.9;%%�������ѡȡ�ĸ���?
    % NI = 30000 * DIM / 10;%%������
    BW = 0.01;%%΢����С
    PAR = 0.3;%%΢���ĸ���
    ansr = zeros(HMS , DIM);%%����⣬HMS*DIM��С�ľ���
    F = zeros(HMS , 1);%%��¼Ŀ�꺯����?
    FP = zeros(1 , NI);%%��¼������м���������������Ŀ�꺯���?
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
    range = [Xmin, Xmax];%%��ֵ��ȡֵ��Χ

    
    %%��ʼ�������?
    for i = 1 : HMS
        for j = 1 : DIM
            ansr(i, j) = getrand(range(1) , range(2));%%����ʼ��ÿ������
        end
    end
    [fitval, fitindex] = min(F);
    F = objfunc.obj_value(ansr');%%��¼ÿ��������Ŀ�꺯����?
    globalBestValue = fitval(1,1);   
    globalBestSolution = ansr(fitindex(1),:);
    %% *****************==- Main loop -==**********************
    for T = 1 : NI
        ans_temp = zeros(1, DIM);%%��ʱ��¼���������?

        %%����µ�����?
        for i = 1 : DIM
            pos = rand();
            if (pos < HMCR)
                ans_temp(i) = ansr(getnum(HMS) , i);%%ֱ���ڼ�����ж�Ӧ�������ѡȡ
                posr = rand();
                if (posr < PAR)
                    ans_temp(i) = ans_temp(i) + dt(BW);%%��ѡȡ��ֵ����΢��
                end
            else
                ans_temp(i) = getrand(range(1), range(2));%%�ڿ���������
            end
        end

        [F_worst , num_worst] = max(abs(F));%%�ҳ�Ŀǰ������е�����?num_worst��¼��������

        %%�жϸ��¼����?
        F_temp = objfunc.obj_value(ans_temp');%%�µ�������Ŀ�꺯��ֵ
        if (abs(F_temp) < F_worst)
            F(num_worst) = F_temp;%%Ŀ�꺯��ֵ����
            ansr(num_worst , :) = ans_temp;%%��������
        end
		if globalBestValue >= F_temp
            globalBestValue = F_temp;
            globalBestSolution(1,:) = ans_temp;
        end
        if T > NI
            break;
        end
        
        [F_best , num_best] = min(abs(F));%%�ҳ�Ŀǰ������е�����?num_best��¼��������
        FP(T) = F(num_best); %%��¼������м���������������Ŀ�꺯���?
%         F_best
        if mod(T,100) == 0
            outcome = [ outcome ;[T FP]];
        end
    end


    %% ****************==- �㷨�����ݴ��� -==*********************
    best_fvalue = F_best;
    best_solution = globalBestSolution;
    run_series = outcome;
    Altime = cputime - mytime;                  % ����ʱ��
    run_info = [Altime];
end

    