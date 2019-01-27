classdef f9 < func_class 
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
            f1 = sum(x.^2)/30;
            f1 = -20*exp(-0.2*sqrt(f1));
            f2 = sum(cos(2*pi*x))/length(x);
            f2 = exp(f2);
            f = f1 - f2 + 20 + exp(1);
        end
        
       %% 构造函数 constructed function
        function obj = f9(specific)
            % specific 预留值,使用结构数组,用于传递额外的参数
            % 设置函数参数的默认值
                if ~exist('specific','var')
                    specific = [];  % 参数初始化
                end
                if ~isfield(specific,'Xdim') specific.Xdim = 30;  end % 变量维数
                if ~isfield(specific,'Xmin') specific.Xmin = -5.12*ones(specific.Xdim,1); end  % 下界 
                if ~isfield(specific,'Xmax') specific.Xmax = 5.12*ones(specific.Xdim,1);  end   % 上界
                if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % 目标函数最小值
                if ~isfield(specific,'Tex')  specific.Tex  =  '-20\exp\left[-0.2\sqrt{\frac{1}{D}\sum\limits_{i=1}^Dx_i^2}\right]-\exp\left(\frac{1}{D}\sum_{i=1}^D\cos(2\pix_i)\right)';   end   % 函数tex表达式
                if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % 目标函数值个数
                if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % 不等式约束
                if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % 等式约束
            obj = obj@func_class(specific);
        end
   end
end    
