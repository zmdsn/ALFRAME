classdef unimodal_f5  < func_class 
   methods
       %% 主运算函数
        function  y = fun_fitness(obj, x , specific )
			swarmsize = size( x , 2 );
            if swarmsize > 1
            % 计算一个种群的函数值
				y = [];  
                for i = 1:swarmsize
                    y(i) = fun_fitness(obj, x(:,i) );
                end
            else
				% 计算一个变量的函数值
                n = length(x);
                y = 0;
                for i=1:(n-1)
                    y = y + 100*(x(i+1)-x(i)^2)^2 + (x(i)-1)^2;
                end
            end
        end
       %% 构造函数
        function obj = unimodal_f5(specific)
            % 定义上下界以及维数
            specific.Xdim = 30 ;         %变量维数
            specific.Xmin = -100*ones(specific.Xdim,1);   %下界
            specific.Xmax = 100*ones(specific.Xdim,1);  %上界
            obj = obj@func_class(specific);
        end
   end
end