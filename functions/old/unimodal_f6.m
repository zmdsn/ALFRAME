classdef unimodal_f6  < func_class 
   methods
       %% 主运算函数
        function  y = fun_fitness(obj, x , specific )
            % 注意这里x为列变量
            % specific 预留值,建议使用结构数组,用于传递额外的参数
            swarmsize = size( x , 2 );
            if swarmsize > 1
            % 计算一个种群的函数值
				y = zeros(swarmsize,1);  
                for i = 1:swarmsize
                    y(i) = fun_fitness(obj, x(:,i) );
                end
            else
            % 计算一个变量的函数值
				ff = floor(x+0.5).^2;
                y = sum(ff);  
            end
        end
       %% 构造函数
        function obj = unimodal_f6(specific)
            % 定义上下界以及维数
            specific.Xdim = 30;  %变量维数
            specific.Xmin = -100*ones(specific.Xdim,1); %下界
            specific.Xmax = 100*ones(specific.Xdim,1);  %上界
            obj = obj@func_class(specific);
        end
   end
end