classdef mof1 < mobjfunc_class 
   methods
       %% 主运算函数
        function  y = fun_fitness(obj, x , specific )
            % 注意这里x为列变量
            % specific 预留值,建议使用结构数组,用于传递额外的参数
            swarmsize = size( x , 2 );
            if swarmsize > 1
            % 计算一个种群的函数值
				y = zeros(swarmsize,obj.objNum);  
                for i = 1:swarmsize
                    y(i,:) = fun_fitness(obj, x(:,i) );
                end
            else
            % 计算一个变量的函数值
                y = [];
                y(1) = x(1)^2;
                y(2) = (x(1)+2)^2;
            end
            if size(y,2) ~= obj.objNum
                warning('The number of decision variables does not match . Kindly check your objective function')
            end
        end
        
       %% 构造函数
        function obj = mof1(in_spec)
            % 构造函数 ,构造一个函数类
            dimension = 3;  % 变量维数
            minimum = -4*ones(dimension,1); % 下界 
            maximum = 4*ones(dimension,1);  % 上界
            specific.Fmin = 'x\in [0 2]';
            obj_num = 2;
            specific.Tex = 'f(1)=x^2\\f(2)=(x-2)^2';  % 函数tex表达式
            switch nargin
                case 0
                    if length(obj_num) == 1 && obj_num > 1
                        specific.objNum = obj_num;
                    else
                        error('Please enter the correct parameters : objNum');
                    end
                    if length(dimension) == 1
                        specific.Xdim = dimension;
                    else
                        error('Please enter the correct parameters : Xdim');
                    end
                    if length(maximum) == 1 
                        specific.Xmax = maximum*ones(dimension,1);    % 上界    
                    else
                        specific.Xmax = maximum;    % 上界                        
                    end
                    if length(minimum) == 1 
                        specific.Xmin = minimum*ones(dimension,1);    % 下界    
                    else
                        specific.Xmin = minimum;    % 下界                        
                    end    
                case 1
                    if isfield(in_spec,'Xdim') 
                        if length(in_spec.Xdim) == 1
                            specific.Xdim = in_spec.Xdim;
                        else
                            error('Please enter the correct parameters : Xdim');
                        end
                    elseif exist('dimension')
                        specific.Xdim = dimension;
                    end
                    if isfield(in_spec,'Xmin')
                        if length(in_spec.Xmin) == 1
                            specific.Xmin = in_spec.Xmin*ones(in_spec.Xdim,1);
                        else
                            specific.Xmin = in_spec.Xmin;
                        end
                    else
                        if length(maximum) == 1 
                            specific.Xmin = maximum*ones(dimension,1);    % 上界    
                        else
                            specific.Xmax = maximum;    % 上界                        
                        end
                    end
                    if isfield(in_spec,'Xmax') 
                        if length(in_spec.Xmax) == 1
                            specific.Xmax = in_spec.Xmax*ones(in_spec.Xdim,1);
                        else
                            specific.Xmax = in_spec.Xmax;
                        end
                    else
                        if length(minimum) == 1 
                            specific.Xmin = minimum*ones(dimension,1);    % 下界    
                        else
                            specific.Xmin = minimum*ones(dimension,1);    % 下界                        
                        end                            
                    end
                    if isfield(in_spec,'objNum') 
                        if length(in_spec.objNum) == 1 && in_spec.objNum >= 1
                            specific.objNum = in_spec.objNum;
                        else
                            error('Please enter the correct parameters : objNum');
                        end
                    elseif exist('obj_num')
                        specific.objNum = obj_num;
                    end
            end
            obj = obj@mobjfunc_class(specific);
        end
   end
end