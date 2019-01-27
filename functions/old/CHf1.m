classdef CHf1 < func_class 
   methods
        %% �����㺯��
        function [f, g, h] = object_function(obj, x , specific )
            % ע������xΪ������
            % specific Ԥ��ֵ,����ʹ�ýṹ����,���ڴ��ݶ���Ĳ���?            % f(x) = 
            % g(x) <= 0
            % h(x) = 0
            g = zeros(obj.Fneq,1);
            h = zeros(obj.Feq,1);
			f = sum((x-50).^2) ;
        end
        
       %% ���캯��
        function obj = CHf1(specific)
            % specific Ԥ��ֵ,ʹ�ýṹ����,���ڴ��ݶ���Ĳ���?            % ���ú�������Ĭ��ֵ
            if ~exist('specific','var')
                specific = [];  % �����ʼ��
            end
            if ~isfield(specific,'Xdim') specific.Xdim = 30;  end % ����ά��
            if ~isfield(specific,'Xmin') specific.Xmin = -100*ones(specific.Xdim,1); end  % �½� 
            if ~isfield(specific,'Xmax') specific.Xmax = 100*ones(specific.Xdim,1);  end   % �Ͻ�
            if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % Ŀ�꺯����Сֵ
            if ~isfield(specific,'Tex')  specific.Tex  = 'f=\sum\limits_{i=1}^D x_i^2';  end   % ����tex����?            if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % Ŀ�꺯��ֵ����
            if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % ����ʽԼ��g
            if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % ��ʽԼ��h
            obj = obj@func_class(specific);
        end
   end
end