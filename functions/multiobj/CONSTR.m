classdef  CONSTR < func_class 
   methods
       %% 主运算函数
        function [f, g, h] = object_function(obj, x , specific )
            % 注意这里x为列向量
            % specific 预留值,建议使用结构数组,用于传递额外的参数
            % f(x) = 
            % g(x) <= 0
            % h(x) = 0
            % CEC07 输入为行向量,输出为行向量,注意转换
            g = zeros(obj.Fneq,1);
            h = zeros(obj.Feq,1);
            f = [0;0];
			f(1) = x(1);
            f(2) = (1+x(2))/x(1);
            g(1) = 6 - x(2) - 9*x(1);
            g(2) = 1 + x(2) - 9*x(1);
        end
        
       %% 构造函数 constructed function
        function obj = CONSTR(specific)
            % specific 预留值,使用结构数组,用于传递额外的参数
            % 设置函数参数的默认值
            if ~exist('specific','var')
                specific = [];  % 参数初始化
            end
            if ~isfield(specific,'Xdim') specific.Xdim = 2;  end % 变量维数
            if ~isfield(specific,'Xmin') specific.Xmin = [0.1 0]; end  % 下界 
            if ~isfield(specific,'Xmax') specific.Xmax = [1 5];  end   % 上界
            if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % 目标函数最小值
            if ~isfield(specific,'Tex')  specific.Tex  = 'f_1(x)=x^2\\f_2(x)=(x-2)^2';  end   % 函数tex表达式
            if ~isfield(specific,'Nobj') specific.Nobj = 2;  end  % 目标函数值个数
            if ~isfield(specific,'Fneq') specific.Fneq = 2;  end  % 不等式约束
            if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % 等式约束
            obj = obj@func_class(specific);
        end
   end
end