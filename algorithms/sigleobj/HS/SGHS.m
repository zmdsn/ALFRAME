function [best_fvalue,best_solution,run_series,run_info,user_set] = SGHS(func_name , algo_config,func_config,user_config)
	getrand = @(l, r) l + (r - l) .* rand();%%�ڼ������ȡ
	getnum = @(num) fix(1 + rand() * num);%%�ڼ������ȡ
	dt = @(num) 2 * num * (rand() - 0.5);%%�������΢��
    mytime = cputime;
	global BW_MAX BW_MIN N
    
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
% 	N = objfunc.Xdim;%% ά��
%     Xmin =  objfunc.Xmin(1);
%     Xmax = objfunc.Xmax(1);
	range = [Xmin, 10];%%��ֵ��ȡֵ��Χ
	HMCRm = 0.98;%%�������ѡȡ�ĸ��ʵĳ�ʼ��ֵ
	HMCRsd = 0.001;%%�������ѡȡ�ĸ��ʵı�׼���С
	PARm = 0.9;%%΢���ĸ��ʵĳ�ʼ��ֵ
	PARsd = 0.01;%%΢���ĸ��ʵı�׼���С
	LP = 100;%%ѧϰ���ڳ���

	ansr = zeros(HMS, N);%%����⣬HMS*N��С�ľ���
	
	F = zeros(1, HMS);%%��¼Ŀ�꺯���С

	BW_MAX = (range(2) - range(1)) / 20;%%΢����С
	BW_MIN = 0.0005;

	%%��ʼ�������
	for i = 1 : HMS
		for j = 1 : N
			ansr(i, j) = getrand(range(1), range(2));%%����ʼ��ÿ������
        end
		F(i) = [objfunc.obj_value( ansr(i, :)')]'; %%��¼ÿ��������Ŀ�꺯���ֵ
    end
    [fitval, fitindex] = min(F);
    globalBestValue = fitval(1,1);   
    globalBestSolution = ansr(fitindex(1),:);
    outcome = [HMS globalBestValue];
    
	ANSP = zeros(NI, N);%%��¼������в����ÿ������
	FP = zeros(1, NI);%%��¼������в����ÿ��������Ŀ�꺯��ֵ
	HMCR_temp = zeros(1, LP);%%��¼��ѡ�е������Ķ�Ӧ��HMCR��ֵ
	PAR_temp = zeros(1, LP);%%��¼��ѡ�е������Ķ�Ӧ��PAR��ֵ
	t = 1;%%��¼��¼��ѡ�е������ĸ���ֵ
	lp = 1;%%��ʼ��ѧϰ���� 

	%%���ѭ��
	for T = 1 : NI
		ans_temp = zeros(1, N);%%��ʱ��¼���������
		
		%%����ѧϰ�Ͳ���HMCR��PAR��ֵ
		if(lp > LP)
			if(t > 1)   
				HMCRm = mean(HMCR_temp(1:t-1));%%��ѡ�е�ȡ��ֵ
				PARm = mean(PAR_temp(1:t-1));%%��ѡ�е�ȡ��ֵ
			end
			lp = 1;%%���³�ʼ��
			t = 1;%%���³�ʼ��
		end
		HMCR =  normrnd(HMCRm , HMCRsd);
		PAR =  normrnd(PARm , PARsd);
		
		%%�ҳ�Ŀǰ������е���������Ž�
		num_i1 = 1;%%num_i1��¼��������
		num_i2 = 1;%%num_i2��¼���Ž�����
		for i = 1 : HMS
			if (F(num_i1) < F(i))
				num_i1 = i;
			end
			if (F(num_i2) > F(i))
				num_i2 = i;
			end
		end
		
		%%����µ�����
		for i = 1 : N
			pos = rand();
			if (pos < HMCR)
				ans_temp(i) = ansr(getnum(HMS), i) + dt(BW(T, NI));%%ֱ���ڼ�����ж�Ӧ�������ѡȡ������΢��
				posr = rand();
				if (posr < PAR)
					ans_temp(i) = ansr(num_i2 , i);%%ֱ��ѡȡ���ŵ�һ����Ӧ��ֵ
				end
			else
				ans_temp(i) = getrand(range(1), range(2));%%�ڿ���������
			end
		end

		%%�жϸ��¼����
		F_temp = [objfunc.obj_value(ans_temp')]';%%�µ�������Ŀ�꺯��ֵ
		if (F_temp < F(num_i1))
			F(num_i1) = F_temp;%%Ŀ�꺯��ֵ����
			ansr(num_i1, :) = ans_temp;%%��������
			HMCR_temp(t) = HMCR;%%��¼��Ӧ��HMCR��ֵ
			PAR_temp (t) = PAR;%%��¼��Ӧ��PAR��ֵ
			t = t + 1;
		end
		if (globalBestValue >= F_temp)
            globalBestValue = F_temp;
            globalBestSolution(1,:) = ans_temp;
        end
        if(mod(T,100) == 0)
            outcome = [outcome ; (T + HMS) globalBestValue];
        end 
		lp = lp + 1;%%ѧϰ�����1
		%%�洢�������ɵ����� 
		ANSP(T, :) = ans_temp;
		FP(T) = F_temp;
	end

	% plot3(ANSP(1 : 10 : NI, 1), ANSP(1 : 10 : NI, 2), ANSP(1 : 10 : NI, 3), 'r.', 'markersize', 1);
% 
% 	ansr(:, 1)'
% 	ansr(:, 2)'
% 	ansr(:, 3)'
    %% ****************==- �㷨�����ݴ��� -==*********************
    best_fvalue = globalBestValue;
    best_solution = globalBestSolution;
    run_series = outcome;
    Altime = cputime - mytime;                  % ����ʱ��
    run_info = [Altime];
    
    
%     altime = cputime - mytime;
%     display([ '       SGHS  ' num2str(altime)]);
% 
% %     display('      SGHS');
%     display(['      ����ֵ��' num2str(globalBestValue)]);
end

function result = BW(T, NI)

	global BW_MAX BW_MIN

	tmp = NI / 2;

	if (T < tmp)
		result = BW_MAX - 2 * T * (BW_MAX - BW_MIN) / NI;
	else
		result = BW_MIN;
	end
end
    