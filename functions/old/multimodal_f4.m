classdef multimodal_f4 < func_class 
   methods
       %% 主运算函数
        function  y = fun_fitness(obj, x , specific )
            % 注意这里x为列变量
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
              y = 0;  
              n = length(x);
              f1 = sum(x.^2)/4000;
              f2 = prod(cos(x'./sqrt(1:n)));
              y = f1 - f2 + 1;
            end
        end
       %% 构造函数
        function obj = multimodal_f4(specific)
            % 定义上下界以及维数
            specific.Xdim = 30;  %变量维数
            specific.Xmin = -32*ones(specific.Xdim,1); %下界
            specific.Xmax = 32*ones(specific.Xdim,1);  %上界
            obj = obj@func_class(specific);
        end
   end
end