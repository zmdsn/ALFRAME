function [best_fvalue,best_solution,run_series,run_info,user_set] = OLPSO(func_name , algo_config,func_config,user_config)

%**************************************************************************************************
%Reference:  Zhan, Zhi Hui, J. Zhang, and O. Liu. "Orthogonal learning particle swarm optimization." 
%              Conference on Genetic and Evolutionary Computation ACM, 2011:1763-1764.
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
    default_set = struct('PopuSize' , 50,...
                         'maxFES'   , 300000  );
                     
    if ~exist('func_name','var') func_name = @sin; end
    if ~exist('algo_config','var') algo_config = []; end
    if ~exist('func_config','var') func_config = []; end
    % get parameters
    [algo_config, func_config, objfunc, PopuSize, maxFES, DIM, Xmin, Xmax, f, init] = ...
        get_parameters( func_name , algo_config , func_config,default_set );
    mytime = cputime;
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);
        
    %% ****************==- main body  -==***********************
  
        dim = DIM;
        SS = PopuSize;
        w = 0.725;  c  =2.0;
        vmax = (Xmax-Xmin)/5;
        FES = 0;  outcome = [];
        it = zeros(SS,1);
        swarm=rand(dim,SS)*(Xmax-Xmin)+Xmin;
        vel=rand(dim,SS)*2*vmax-vmax;
        bestpos=swarm;    sswarm=swarm;
        fswarm = zeros(1,SS);  fbestpos = fswarm;
        vvel=vel;   
        for i=1:SS
            fswarm(i) = f(swarm(:,i),func_config);
            
            fbestpos(i) = fswarm(i); 
            FES = FES + 1;
        end
        [fxopt,g_over] = min(fbestpos);
        xopt = bestpos(:, g_over);  
        outcome = [outcome; [FES fxopt]];      
        % construct orthogonal array
        m=2^ceil(log(dim+1)/log(2));
        L=ones(m,m-1);
        u=log(m)/log(2);
        for k=1:u
              bb=2^(k-1);
              for aa=1:m
                  L(aa,bb)=mod(floor((aa-1)/2^(u-k)),2);
              end
        end
        for i=2:u
           bb=2.^(0:(u-1));
           s=nchoosek(bb,i);
           ss=sum(s,2);
           for j=1:size(s,1)
             LL=zeros(m,1);
              for k=1:size(s,2)
                  LL=LL+L(:,s(j,k));
              end
             L(:,ss(j))=mod(LL,2); 
           end
        end
        L=L(:,randperm(dim));
           
           % orthogonal learning        
         for i=1:SS 
             for j=1:m
                 for k=1:dim
                     if(L(j,k)==0)
                         p0(k)=xopt(k);
                     else
                         p0(k)=bestpos(k,i);
                     end
                 end
                 F(j) = f(p0',func_config);
                 FES = FES + 1;
             end
             for ii=1:dim
                 l=L(:,ii);
                 [oo,p]=min([sum(F)-F*l,F*l]);
                 ll(ii)=p-1;
             end
             for k=1:dim
                 if(ll(k)==0)
                     p0(k)=xopt(k);
                 else
                     p0(k)=bestpos(k,i);
                 end
             end
             PO(:,i)=p0';
         end
           
        while FES<maxFES
            for i=1:SS
                R2=rand(dim,1);
                vel(:,i)=w.*vel(:,i)+c.*R2.*(PO(:,i)-swarm(:,i));
            end
            vel(vel>vmax) = vvel(vel>vmax);
            vel(vel<-vmax) = vvel(vel<-vmax);
            swarm = swarm + vel;   
            vvel = vel;
            swarm(swarm>Xmax) = sswarm(swarm>Xmax);
            swarm(swarm<Xmin) = sswarm(swarm<Xmin);
            sswarm = swarm;

           for i=1:SS
               fswarm(i) = f(swarm(:,i),func_config); %benchmark_func(swarm(:,i)', func_num, o, A, M, a, alpha, b);
               FES = FES + 1;
           end
           bestpos(:,fswarm<fbestpos) = swarm(:,fswarm<fbestpos);
           fbestpos(fswarm<fbestpos) = fswarm(fswarm<fbestpos);
           [fxopt,g_over] = min(fbestpos);
           xopt = bestpos(:,g_over);

            for i=1:SS
                if(fswarm(i)>fxopt)
                    it(i)=it(i)+1;
                end
            end

            for i=1:SS
                if(it(i)>5)
                    it(i)=0;
                    for j=1:m
                         % orthogonal learning
                         for k=1:dim
                             if(L(j,k)==0)
                                 p0(k)=xopt(k);
                             else
                                 p0(k)=bestpos(k,i);
                             end
                         end
                         F(j) = f(p0',func_config); %benchmark_func(p0, func_num, o, A, M, a, alpha, b);
                         FES = FES + 1;
                     end
                     for ii=1:dim
                         l=L(:,ii);
                         [oo,p]=min([sum(F)-F*l,F*l]);
                         ll(ii)=p-1;
                     end
                     for k=1:dim
                         if(ll(k)==0)
                             p0(k)=xopt(k);
                         else
                             p0(k)=bestpos(k,i);
                         end
                     end
                     PO(:,i)=p0';
                end
            end 
            if FES >= maxFES
                break;          % to maintain data length consistent
            else
                outcome = [outcome; [FES fxopt]];                
            end
        end  
 

    %% ****************==- collating the results -==*********************
    run_series = outcome;
    Altime = cputime - mytime;                  
    run_info = [Altime];
    best_fvalue = fxopt;
    best_solution = xopt;
end
