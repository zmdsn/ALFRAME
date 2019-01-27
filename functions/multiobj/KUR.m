classdef  KUR < func_class 
%**************************************************************************
% References
    % [1] *Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal, and T. Meyarivan*, |A Fast
    % Elitist Multiobjective Genetic Algorithm: NSGA-II|, IEEE Transactions on 
    % Evolutionary Computation 6 (2002), no. 2, 182 ~ 197.
    % [2] *N. Srinivas and Kalyanmoy Deb*, |Multiobjective Optimization Using 
    % Nondominated Sorting in Genetic Algorithms|, Evolutionary Computation 2 
    % (1994), no. 3, 221 ~ 248.
% Author  : unclear
% Adapter : Algori
% Email : zmdsn@126.com
% programmed: Sept 29, 2016
%**************************************************************************
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
            x1 = x(1:end-1);
            x2 = x(2:end);
			f(1) = sum(-10*exp(-0.2.*sqrt(x1.^2+x2.^2)));
            f(2) = sum(abs(x).^0.8+5*sin(x.^3));
        end
        
       %% 构造函数 constructed function
        function obj = KUR(specific)
            % specific 预留值,使用结构数组,用于传递额外的参数
            % 设置函数参数的默认值
            if ~exist('specific','var')
                specific = [];  % 参数初始化
            end
            if ~isfield(specific,'Xdim') specific.Xdim = 3;  end % 变量维数
            if ~isfield(specific,'Xmin') specific.Xmin = -5*ones(specific.Xdim,1); end  % 下界 
            if ~isfield(specific,'Xmax') specific.Xmax = 5*ones(specific.Xdim,1);  end   % 上界
            if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % 目标函数最小值
            if ~isfield(specific,'Tex')  specific.Tex  = 'f_1(x)=x^2\\f_2(x)=(x-2)^2';  end   % 函数tex表达式
            if ~isfield(specific,'Nobj') specific.Nobj = 2;  end  % 目标函数值个数
            if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % 不等式约束
            if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % 等式约束
            obj = obj@func_class(specific);
        end
   end
end