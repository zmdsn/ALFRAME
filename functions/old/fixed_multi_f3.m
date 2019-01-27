classdef fixed_multi_f3 < func_class 
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
       %% 主运算函数
        function  [y,g,h] = object_function(obj, x , specific )
            % 注意这里x为列变量
            g = 0;
            h = 0;
            n = length(x);
            F = 1 + (x+1)*0.25;
            y = 10*sin(pi*F(1))^2 + (F(n)-1)^2;
            for i=1:(n-1)
                y = y + (F(i)-1)^2*(1+10*sin(pi*F(i+1))^2);
            end
            y = pi*y/n;
            for i=1:n
                y = y + obj.u(x(i),10,100,4);
            end
        end
       %% 构造函数
        function obj = fixed_multi_f3(specific)
            % 定义上下界以及维数
            % 非常有趣的一个函数
            if ~exist('specific','var')
                specific = [];  % 参数初始化
            end
            if ~isfield(specific,'Xdim') specific.Xdim = 2;  end % 变量维数
            if ~isfield(specific,'Xmin') specific.Xmin = -10*ones(specific.Xdim,1); end  % 下界 
            if ~isfield(specific,'Xmax') specific.Xmax = 10*ones(specific.Xdim,1);  end   % 上界
            if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % 目标函数最小值
            if ~isfield(specific,'Tex')  specific.Tex  =  'f=\frac{1}{4000}\sum\limits_{i=1}^Dx_i^2-\prod_{i=1}^D\cos\left(\frac{x_i}{\sqrt{i}}\right)';   end   % 函数tex表达式
            if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % 目标函数值个数
            if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % 不等式约束
            if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % 等式约束
            obj = obj@func_class(specific);
        end
   end
end