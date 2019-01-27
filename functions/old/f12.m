classdef f12  < func_class 
   methods
        function FF = ux(~,x,a,k,m)
			FF = 0;
			if (x>a) FF = k*power(x-a,m);  end
			if (x<-a) FF = k*power(-x-a,m);  end
        end
       %% 主运算函数
%         function  y = fun_fitness(obj, x , specific )
        function [f, g, h] = object_function(obj, x , specific )
            % 注意这里x为列变量
			swarmsize = size( x , 2 );
            if swarmsize > 1
            % 计算一个种群的函数值
				y = [];  
                for i = 1:swarmsize
                    y(i) = fun_fitness(obj, x(:,i) );
                end
            else
				% 计算一个变量的函数值
                y = 0;
                n = length(x);
                y = sin(3*pi*x(1))^2 + (x(n)-1)^2*(1+sin(2*pi*x(n))^2);
                for i=1:(n-1)
                    y = y + (x(i)-1)^2*(1+sin(3*pi*x(i+1))^2);
                end
                y = 0.1*y;
                for i=1:n
                    y = y + ux(obj,x(i),5,100,4);
                end
            end
            f = y;
            g = [];
            h = [];
        end
       %% 构造函数
%         function obj = f12(specific)
%             % 定义上下界以及维数
%             specific.Xdim = 30;  %变量维数
%             specific.Xmin = -50*ones(specific.Xdim,1); %下界
%             specific.Xmax = 50*ones(specific.Xdim,1);  %上界
%             obj = obj@func_class(specific);
%         end
        function obj = f12(specific)
            % specific 预留值,使用结构数组,用于传递额外的参数
            % 设置函数参数的默认值
                if ~exist('specific','var')
                    specific = [];  % 参数初始化
                end
                if ~isfield(specific,'Xdim') specific.Xdim = 30;  end % 变量维数
                if ~isfield(specific,'Xmin') specific.Xmin = -50*ones(specific.Xdim,1); end  % 下界 
                if ~isfield(specific,'Xmax') specific.Xmax = 50*ones(specific.Xdim,1);  end   % 上界
                if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % 目标函数最小值
                if ~isfield(specific,'Tex')  specific.Tex  =  'f=\frac{1}{4000}\sum\limits_{i=1}^Dx_i^2-\prod_{i=1}^D\cos\left(\frac{x_i}{\sqrt{i}}\right)';   end   % 函数tex表达式
                if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % 目标函数值个数
                if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % 不等式约束
                if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % 等式约束
            obj = obj@func_class(specific);
        end
   end
end