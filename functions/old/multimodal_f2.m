classdef multimodal_f2  < func_class 
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
                y = 0;
                y = sum( x.^2 - 10*cos(2*pi*x) + 10 );
            end
        end
       %% 构造函数
        function obj = multimodal_f2(specific)
            % 定义上下界以及维数
            specific.Xdim = 30  ;        %变量维数
            specific.Xmin = -5.12*ones(specific.Xdim,1);   %下界
            specific.Xmax = 5.12*ones(specific.Xdim,1);  %上界
            obj = obj@func_class(specific);
        end
   end
end