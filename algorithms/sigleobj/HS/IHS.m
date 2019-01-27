function [best_fvalue,best_solution,run_series,run_info,user_set] = IHS(func_name , algo_config,func_config,user_config)
    % fans = 0;%%�Ż�Ŀ��
    % f = @(x, y, z) abs((1 - x) .^ 2 + 100 * (y - x .^ 2) .^ 2 + (y - z .^ 3) .^ 2-fans);%%�Ż���Ŀ�꺯��
    getrand = @(l, r) l + (r - l) .* rand();%%�ڼ������ȡ
    getnum = @(num) fix(1 + rand() * num);%%�ڼ������ȡ
    dt = @(num) 2 * num * (rand() - 0.5);%%�������΢��
    mytime = cputime;
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
    HMCR = 0.9;%%�������ѡȡ�ĸ���

    ansr = zeros(HMS, N);%%����⣬HMS*N��С�ľ���
	range = [Xmin, Xmax];%%��ֵ��ȡֵ��Χ
    F = zeros(1, HMS);%%��¼Ŀ�꺯���С

    global BW_MAX BW_MIN PAR_MAX PAR_MIN
    BW_MAX = (range(2) - range(1)) / 20;%%΢����С
    BW_MIN = 0.0001;
    PAR_MAX = 0.99;%%΢���ĸ���
    PAR_MIN = 0.01;

    %%��ʼ�������
    for i = 1 : HMS
        for j = 1 : N
            ansr(i, j) = getrand(range(1), range(2));%%����ʼ��ÿ������
        end
        F(i) = objfunc.obj_value( ansr(i, :)');%%��¼ÿ��������Ŀ�꺯���ֵ
    end
    [fitval, fitindex] = min(F);
    globalBestValue = fitval(1,1);   
    globalBestSolution = ansr(fitindex(1),:);
    outcome = [HMS globalBestValue];

    ANSP = zeros(NI, N);%%��¼������в����ÿ������
    FP = zeros(1, NI);%%��¼������в����ÿ��������Ŀ�꺯��ֵ

    %%���ѭ��
    for T = 1 : NI
        ans_temp = zeros(1, N);%%��ʱ��¼���������

        %%����µ�����
        for i = 1 : N
            pos = rand();
            if (pos < HMCR)
                ans_temp(i) = ansr(getnum(HMS), i);%%ֱ���ڼ�����ж�Ӧ�������ѡȡ
                posr = rand();
                if (posr < PAR(T, NI))
                    ans_temp(i) = ans_temp(i) + dt(BW(T, NI));%%��ѡȡ��ֵ����΢��
                end
            else
                ans_temp(i) = getrand(range(1), range(2));%%�ڿ���������
            end
        end

        %%�ҳ�Ŀǰ������е�����
        num_i = 1;%%num_i��¼��������
        for i = 1 : HMS
            if (F(num_i) < F(i))
                num_i = i;
            end
        end

        %%�жϸ��¼����
        F_temp = [objfunc.obj_value( ans_temp')]';%%�µ�������Ŀ�꺯��ֵ
        if (F_temp < F(num_i))
            F(num_i) = F_temp;%%Ŀ�꺯��ֵ����
            ansr(num_i, :) = ans_temp;%%��������
        end
		if (globalBestValue >= F_temp)
            globalBestValue = F_temp;
            globalBestSolution(1,:) = ans_temp;
        end
        if(mod(T,100) == 0)
            outcome = [outcome ; (T + HMS) globalBestValue];
        end 
       %%�洢�������ɵ����� 
       ANSP(T, :) = ans_temp;
       FP(T) = F_temp;
    end

    %% ****************==- �㷨�����ݴ��� -==*********************
    best_fvalue = globalBestValue;
    best_solution = globalBestSolution;
    run_series = outcome;
    Altime = cputime - mytime;                  % ����ʱ��
    run_info = [Altime];
end

function result = PAR(T, NI)
    global PAR_MAX PAR_MIN
    result = PAR_MIN + T * (PAR_MAX - PAR_MIN) / NI;
end

function result = BW(T, NI)
    global BW_MAX BW_MIN

    % ans = BW_MIN + (BW_MAX - BW_MIN) * exp(1-T);
    result = BW_MAX * exp((log(BW_MIN/BW_MAX)/NI) * T);
end