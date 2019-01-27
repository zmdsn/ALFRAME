classdef f14  < func_class 
   methods
       %% 主运算函数
        function  [f, g, h] = object_function(obj, x , specific )
			swarmsize = size( x , 2 );
            if swarmsize > 1
            % 计算一个种群的函数值
				y = [];  
                for i = 1:swarmsize
                    y(i) = object_function(obj, x(:,i) );
                end
            else
				% 计算一个变量的函数值
                a=[0.1957,0.1947,0.1735,0.16,0.0844,0.0627,0.0456,0.0342,0.0323,0.0235,0.0246];
                b=1./[0.25,0.5,1,2,4,6,8,10,12,14,16];
                y = 0; n=length(a);
                for i=1:n
                   y = y + (a(i) - (x(1)*b(i)*(b(i) + x(2)))/(b(i)*(b(i)+x(3)) + x(4)))^2;
                end
            end    
            f = y;
            g = [];
            h = [];
        end
       %% 构造函数
%         function obj = f14(specific)
%             % 定义上下界以及维数
%             specific.Xdim = 30;          %变量维数
%             specific.Xmin = -65.536*ones(specific.Xdim,1);   %下界
%             specific.Xmax = 65.536*ones(specific.Xdim,1);  %上界
%             obj = obj@func_class(specific);
%         end
        function obj = f14(specific)
            % specific 预留值,使用结构数组,用于传递额外的参数
            % 设置函数参数的默认值
                if ~exist('specific','var')
                    specific = [];  % 参数初始化
                end
                if ~isfield(specific,'Xdim') specific.Xdim = 30;  end % 变量维数
                if ~isfield(specific,'Xmin') specific.Xmin = -50*ones(specific.Xdim,1); end  % 下界 
                if ~isfield(specific,'Xmax') specific.Xmax = 50*ones(specific.Xdim,1);  end   % 上界
                if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % 目标函数最小值
                if ~isfield(specific,'Tex')  specific.Tex  =  'f=\frac{1}{4000}\sum\limits_{i=1}^Dx_i^2-\prod_{i=1}^D\cos\left(\frac{x_i}{\sqrt{i}}\right)';   end   % 函数tex表达式
                if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % 目标函数值个数
                if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % 不等式约束
                if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % 等式约束
            obj = obj@func_class(specific);
        end
   end
end