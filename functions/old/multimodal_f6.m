classdef multimodal_f6 < func_class 
   methods
        function uout = u(~,x,a,k,m)
            if x > a
                uout = k*(x - a)^m;
            elseif x < -a
                uout = k*(-x - a)^m;
            else
                uout = 0;
            end
        end
        
        function y = value_y(~,x)
            y = 1+(x+1)./4;
        end
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
              y = 0; n = length(x);
              y2 = 0;
              y3 = 0;
              y4 = 0;
              yval = obj.value_y(x);
              y1 = sin(3*pi*x(1))^2;
              for i = 1:n
                y2 = y2 + ((x(i)-1)^2)*( ( 1 + (sin(3*pi*x(i) + 1))^2 ));
                y4 = y4 + obj.u(x(i),10,100,4);
              end
              y3 = ((x(n)-1)^2)*( ( 1 + (sin(2*pi*x(n)))^2 ) );
               
              y = 0.1*( y1 + y2 + y3 ) + y4;

%                  y = 0;
%                 n = length(x);
%                 y = sin(3*pi*x(1))^2 + (x(n)-1)^2*(1+sin(2*pi*x(n))^2);
%                 for i=1:(n-1)
%                     y = y + (x(i)-1)^2*(1+sin(3*pi*x(i+1))^2);
%                 end
%                 y = 0.1*y;
%                 for i=1:n
%                     y = y + obj.u(var(i),5,100,4);
%                 end
            end
        end
       %% 构造函数
        function obj = multimodal_f6(specific)
            % 定义上下界以及维数
            specific.Xdim = 30;  %变量维数
            specific.Xmin = -5*ones(specific.Xdim,1); %下界
            specific.Xmax = 5*ones(specific.Xdim,1);  %上界
            obj = obj@func_class(specific);
        end
   end
end