function [best_fvalue,best_solution,run_series,run_info,user_set] = CLPSO(func_name , algo_config,func_config,user_config)

%**************************************************************************
% Abstract :  
% Susilo, Andy. "COMPREHENSIVE LEARNING PARTICLE SWARM OPTIMIZER (CLPSO) 
% FOR GLOBAL OPTIMIZATION OF MULTIMODAL FUNCTIONS." Undergraduate Theses (2008).
%
% Author : Algo
% Email : zmdsn@126.com
% Date : 9/28/2016
%**************************************************************************

    %% ****************==- Initialization settings -==***********************
    format long;
    format compact;
    rng('shuffle'); 
    % set default
    default_set = struct('PopuSize' , 100,...
                         'maxFES'   , 300000  );
                     
    if ~exist('func_name','var') func_name = @sin; end
    if ~exist('algo_config','var') algo_config = []; end
    if ~exist('func_config','var') func_config = []; end
    % get parameters
    [algo_config, func_config, objfunc, ps, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
        
    %% ****************==- main body  -==***********************
    DD = {};
    ee = 1;
        me = maxFES/ps;
    %     ps = 40;
        outcome = [];
        lu = [Xmin * ones(1, DIM); Xmax * ones(1, DIM)];
        t = 0 : 1 / (ps - 1) : 1;
        t = 5 .* t;
        Pc = 0.0 + (0.5 - 0.0) .* (exp(t) - exp(t(1))) ./ (exp(t(ps)) - exp(t(1)));
        m = 0 .* ones(ps, 1);
        iwt = 0.9 - (1 : me) * (0.7 / me);
        cc = [1.49445 1.49445]; %acceleration constants
        mv = 0.2 * (lu(2, :) - lu(1, :));
        VRmin = repmat(lu(1, :), ps, 1);
        VRmax = repmat(lu(2, :), ps, 1);
        Vmin = repmat(-mv, ps, 1);
        Vmax = -Vmin;

        pos = VRmin + (VRmax - VRmin).* rand(ps, DIM);

        e = f(pos',func_config);

        fitcount = ps;
        vel = Vmin + 2 .* Vmax .* rand(ps, DIM); %initialize the velocity of the particles

        pbest = pos;
        pbestval = e;  %initialize the pbest and the pbest's fitness value

        [gbestval, gbestid] = min(pbestval);
        gbest = pbest(gbestid, :); %initialize the gbest and the gbest's fitness value
        gbestrep = repmat(gbest, ps, 1);

        stay_num = zeros(ps, 1);

        ai = zeros(ps, DIM);
        f_pbest = 1 : ps;
        f_pbest = repmat(f_pbest', 1, DIM);
        for k = 1 : ps
            ar = randperm(DIM);
            ai(k, ar(1 : m(k))) = 1;
            fi1 = ceil(ps * rand(1, DIM));
            fi2 = ceil(ps * rand(1, DIM));

            fi = (pbestval(fi1) < pbestval(fi2)) .* fi1 + (pbestval(fi1) >= pbestval(fi2)) .* fi2;
            bi = ceil(rand(1, DIM) - 1 + Pc(k));

            if bi == zeros(1, DIM),
                rc = randperm(DIM);
                bi(rc(1)) = 1;
            end
            f_pbest(k, :) = bi .* fi + (1 - bi) .* f_pbest(k, :);
        end
        pbest_f = zeros(ps,DIM);
        stop_num = 0;
        for i = 2 : me
            valid = [];
            for k = 1 : ps

                if stay_num(k) >= 5

                    stay_num(k) = 0;
                    ai(k, :) = zeros(1, DIM);
                    f_pbest(k, :) = k .* ones(1, DIM);
                    ar = randperm(DIM);
                    ai(k, ar(1 : m(k))) = 1;
                    fi1 = ceil(ps * rand(1, DIM));
                    fi2 = ceil(ps * rand(1, DIM));
                    fi = (pbestval(fi1) < pbestval(fi2)) .* fi1 + (pbestval(fi1) >= pbestval(fi2)) .* fi2;
                    bi = ceil(rand(1, DIM) - 1 + Pc(k));
                    if bi == zeros(1, DIM),
                        rc = randperm(DIM);
                        bi(rc(1)) = 1;
                    end
                    f_pbest(k, :) = bi .* fi + (1 - bi) .* f_pbest(k, :);
                end

                for dimcnt = 1 : DIM
                    k
                    dimcnt
                    f_pbest(k, dimcnt)
                    pbest(f_pbest(k, dimcnt), dimcnt)
                    pbest_f(k, dimcnt)
                    pbest_f(k, dimcnt) = pbest(f_pbest(k, dimcnt), dimcnt);
                end

                temp = pos(k, :);
                temp_ = temp + vel(k, :);

                aa(k, :) = cc(1) .* (1 - ai(k, :)) .* rand(1, DIM) .* (pbest_f(k, :) - pos(k, :)) + cc(2) .* ai(k, :) .* rand(1, DIM) .* (gbestrep(k, :) - pos(k, :));
                vel(k, :) = iwt(i) .* vel(k, :) + aa(k, :);
                vel(k, :) = (vel(k, :) > mv) .* mv + (vel(k, :) <= mv) .* vel(k, :);
                vel(k, :) = (vel(k, :) < (-mv)) .* (-mv) + (vel(k, :) >= (-mv)) .* vel(k, :);
                pos(k, :) = pos(k, :) + vel(k, :);

                if (sum(pos(k, :) > VRmax(k, :)) + sum(pos(k, :) < VRmin(k, :))) == 0;
                    valid = [valid k];
                    fitcount = fitcount + 1;
                end

            end

            if ~isempty(valid)

                e(valid, 1) = f(pos(valid, :)',func_config);

                for k = 1 : length(valid)

                    tmp = (pbestval(valid(k)) <= e(valid(k)));
                    if tmp == 1
                        stay_num(valid(k)) = stay_num(valid(k)) + 1;
                    end
                    temp = repmat(tmp, 1, DIM);
                    pbest(valid(k), :) = temp .* pbest(valid(k), :) + (1 - temp) .* pos(valid(k), :);
                    pbestval(valid(k)) = tmp .* pbestval(valid(k)) + (1 - tmp) .* e(valid(k)); %update the pbest
                    if pbestval(valid(k)) < gbestval
                        gbest = pbest(valid(k), :);
                        gbestval = pbestval(valid(k));
                        gbestrep = repmat(gbest, ps, 1); %update the gbest
                    end
                end
            end

            if fitcount >= maxFES
                break;
            end
            if (i == me) && (fitcount < maxFES)
                i = i-1;
            end
            DD{ee} = pbest';
            ee = ee+1;
         outcome = [outcome; [i*ps min(pbestval)]];
        end

        
    %% ****************==- collating the results -==*********************
%     save('CLPSOf1.mat','DD')
    run_series = outcome;
    Altime = cputime - mytime;                 
    run_info = [Altime];
    best_fvalue = gbestval;
    best_solution = gbest;
