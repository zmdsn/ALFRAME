classdef f10 < func_class 
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
            f = 0;  
            n = length(x);
            f1 = sum(x.^2)/4000;
            f2 = prod(cos(x'./sqrt(1:n)));
            f = f1 - f2 + 1;
        end
        
       %% 构造函数 constructed function
        function obj = f10(specific)
            % specific 预留值,使用结构数组,用于传递额外的参数
            % 设置函数参数的默认值
            if ~exist('specific','var')
                specific = [];  % 参数初始化
            end
            if ~isfield(specific,'Xdim') specific.Xdim = 30;  end % 变量维数
            if ~isfield(specific,'Xmin') specific.Xmin = -32*ones(specific.Xdim,1); end  % 下界 
            if ~isfield(specific,'Xmax') specific.Xmax = 32*ones(specific.Xdim,1);  end   % 上界
            if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % 目标函数最小值
            if ~isfield(specific,'Tex')  specific.Tex  =  'f=\frac{1}{4000}\sum\limits_{i=1}^Dx_i^2-\prod_{i=1}^D\cos\left(\frac{x_i}{\sqrt{i}}\right)';   end   % 函数tex表达式
            if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % 目标函数值个数
            if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % 不等式约束
            if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % 等式约束
            obj = obj@func_class(specific);
        end
   end
end    
