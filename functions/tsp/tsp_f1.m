classdef tsp_f1 < func_class 
   methods
       %% 主运算函数
        function  [y g h] = object_function(obj,x,specific)
            % 注意这里x为列变量
            % specific 预留值,建议使用结构数组,用于传递额外的参数
            y =0;
            g=[];
            h=[];
            n = obj.Xdim;
            for j = 1:(n - 1)
              y = y + specific.D(x(j),x(j + 1));
            end
            y = y + specific.D(x(n),x(1));
        end
        
       %% 计算城市之间的距离
        function  [n,Distance citys] = cal_distance(obj, specific )
            %% 导入数据
            load ./functions/tsp/data/citys_data.mat
            
            %% 计算城市间相互距离
            n = size(citys,1); % claculate the number of nodes
            Distance = zeros(n,n);
            for i = 1:n
                for j = 1:n
                    if i ~= j
                        Distance(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
                    else
                        Distance(i,j) = 1e-4;      
                    end
                end    
            end
        end
        
       %% 构造函数
        function obj = tsp_f1(in_spec)
            % 构造函数 ,构造一个函数类
            if ~exist('specific','var')
                specific = [];  % 参数初始化
            end
            if ~isfield(specific,'Xdim') specific.Xdim = 30;  end % 变量维数
            if ~isfield(specific,'Xmin') specific.Xmin = -100*ones(specific.Xdim,1); end  % 下界 
            if ~isfield(specific,'Xmax') specific.Xmax = 100*ones(specific.Xdim,1);  end   % 上界
            if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % 目标函数最小值
            if ~isfield(specific,'Tex')  specific.Tex  = 'f=\sum\limits_{i=1}^D x_i^2';  end   % 函数tex表达式
            if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % 目标函数值个数
            if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % 不等式约束
            if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % 等式约束
            obj = obj@func_class(specific);
        end
   end
end