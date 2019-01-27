classdef multimodal_f3  < func_class 
   methods
       %% 主运算函数
        function  y = fun_fitness(obj, x , specific )
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
                y = 0;  n = length(x);
                for i=1:n
                    y = y - x(i)sin(sqrt(abs(x(i))));
                end
            end
        end
       %% 构造函数
        function obj = multimodal_f1(specific)
            % 定义上下界以及维数
            specific.Xdim = 30  ;        %变量维数
            specific.Xmin = -500*ones(specific.Xdim,1);   %下界
            specific.Xmax = 500*ones(specific.Xdim,1);  %上界
            obj = obj@func_class(specific);
        end
   end
end