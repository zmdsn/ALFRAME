classdef multimodal_f8 < func_class 
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
              m = 10;
              m = 2*m;
              y1 = sum((sin(x)).^2);
              y2 = exp(-sum(x.^2));
              y3 = exp( -sum( sin(sqrt(abs(x))).^2 ) );
              y = (y1 - y2)*y3;
            end
        end
       %% 构造函数
        function obj = multimodal_f8(specific)
            % 定义上下界以及维数
            % 非常有趣的一个函数
            specific.Xdim = 30;  %变量维数
            specific.Xmin = -20*ones(specific.Xdim,1); %下界
            specific.Xmax = 20*ones(specific.Xdim,1);  %上界
            obj = obj@func_class(specific);
        end
   end
end