classdef func_handle < func_class 
   methods
        %% Ö÷ÔËËãº¯Êý
        function [f, g, h] = object_function(obj, x , specific )
            % ×¢ÒâÕâÀïxÎªÁÐ±äÁ¿
            % specific Ô¤ï¿½ï¿½Öµ,ï¿½ï¿½ï¿½ï¿½Ê¹ï¿½Ã½á¹¹ï¿½ï¿½ï¿½ï¿½,ï¿½ï¿½ï¿½Ú´ï¿½ï¿½Ý¶ï¿½ï¿½ï¿½Ä²ï¿½ï¿½ï¿?
            % g(x) <= 0
            % h(x) = 0
            g = zeros(obj.Fneq,1);
            h = zeros(obj.Feq,1);
			f = obj.object_function(x);
            if size(f,1) > 1
                warning('too many output,we only use the first one');
                f = f(1);
            end
        end
        
       %% constructed function
        function obj = func_handle(specific)
            % specific 
            % please input column vector
            if ~exist('specific','var')
                specific = [];  % ï¿½ï¿½ï¿½ï¿½ï¿½Ê¼ï¿½ï¿?            
            end
            if ~isfield(specific,'Xdim') specific.Xdim = 30;  end % ï¿½ï¿½ï¿½ï¿½Î¬ï¿½ï¿½
            if ~isfield(specific,'Xmin') specific.Xmin = -100*ones(specific.Xdim,1); end  % ï¿½Â½ï¿½ 
            if ~isfield(specific,'Xmax') specific.Xmax = 100*ones(specific.Xdim,1);  end   % ï¿½Ï½ï¿½
            if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % Ä¿ï¿½êº¯ï¿½ï¿½ï¿½ï¿½Ð¡Öµ
            if ~isfield(specific,'Tex')  specific.Tex  = 'This is a function_handle';  end   % ï¿½ï¿½ï¿½ï¿½texï¿½ï¿½ï¿½Ê?            if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % Ä¿ï¿½êº¯ï¿½ï¿½Öµï¿½ï¿½ï¿½ï¿½
            if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % ï¿½ï¿½ï¿½ï¿½Ê½Ô¼ï¿½ï¿½g
            if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % ï¿½ï¿½Ê½Ô¼ï¿½ï¿½h
            obj = obj@func_class(specific);
        end
   end
end