classdef f11 < func_class 
   methods
        %% 主运算函数
        function [f, g, h] = object_function(obj, x , specific )
            % 注意这里x为列向量
            % specific 预留值,建议使用结构数组,用于传递额外的参数
            % f(x) = 
            % g(x) <= 0
            % h(x) = 0
            g = 0;
            h = 0;
            n = length(x);
            F = 1 + (x+1)*0.25;
            f = 10*sin(pi*F(1))^2 + (F(n)-1)^2;
            for i=1:(n-1)
                f = f + (F(i)-1)^2*(1+10*sin(pi*F(i+1))^2);
            end
            f = pi*f/n;
            for i=1:n
                f = f + ux(obj,x(i),10,100,4);
            end
        end
        
       %% 构造函数 constructed function
        function obj = f11(specific)
            % specific 预留值,使用结构数组,用于传递额外的参数
            % 设置函数参数的默认值
                if ~exist('specific','var')
                    specific = [];  % 参数初始化
                end
                if ~isfield(specific,'Xdim') specific.Xdim = 30;  end % 变量维数
                if ~isfield(specific,'Xmin') specific.Xmin = -600*ones(specific.Xdim,1); end  % 下界 
                if ~isfield(specific,'Xmax') specific.Xmax = 600*ones(specific.Xdim,1);  end   % 上界
                if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % 目标函数最小值
                if ~isfield(specific,'Tex')  specific.Tex  =  'f=\frac{1}{4000}\sum\limits_{i=1}^Dx_i^2-\prod_{i=1}^D\cos\left(\frac{x_i}{\sqrt{i}}\right)';   end   % 函数tex表达式
                if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % 目标函数值个数
                if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % 不等式约束
                if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % 等式约束
            obj = obj@func_class(specific);
        end
        function FF = u(x,a,k,m)
			FF = 0;
			if (x>a) FF = k*power(x-a,m);  end
			if (x<-a) FF = k*power(-x-a,m);  end
        end
        function FF = ux(~,x,a,k,m)
			swarmsize = size( x , 2 );
            n = size( x , 1 );
            y = zeros( 1 , n );
            if swarmsize > 1
            % 计算一个种群的函数值
				y = [];  
                for i = 1:swarmsize
                    y(i) = fun_fitness(obj, x(:,i) );
                end
            else
				% 计算一个变量的函数值
				FF = 0;
                if (x>a) FF = k*power(x-a,m);  end
                if (x<-a) FF = k*power(-x-a,m);  end
            end
        end
   end
end    
