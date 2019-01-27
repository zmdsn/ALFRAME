classdef cec15ni_f15  < func_class 
    methods       
        %% �����㺯��
        function  [f, g, h] = object_function(obj, x , specific )
            % ע������xΪ�б���
            leng = length(x);
            allows = [10, 20, 30];
            idx = find(allows == leng);
            if isempty(idx)
               error(['only allowed dimention is: ' num2str(allows(1:end))]); 
            end
            g = zeros(obj.Fneq,1);
            h = zeros(obj.Feq,1);
            f = cec15_nich_func(x,15);
        end
        %% ���캯��
        function obj = cec15ni_f15(specific)
            % specific Ԥ��ֵ,ʹ�ýṹ����,���ڴ��ݶ���Ĳ���
            % ���ú���������Ĭ��ֵ
            if ~exist('specific','var')
                specific = [];  % ������ʼ��
            end
            if ~isfield(specific,'Xdim') specific.Xdim = 30;  end % ����ά��
            if ~isfield(specific,'Xmin') specific.Xmin = -100*ones(specific.Xdim,1); end  % �½� 
            if ~isfield(specific,'Xmax') specific.Xmax = 100*ones(specific.Xdim,1);  end   % �Ͻ�
            if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % Ŀ�꺯����Сֵ
            if ~isfield(specific,'Tex')  specific.Tex  = 'Function expressions refer to reference.';  end   % ����tex����ʽ
            if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % Ŀ�꺯��ֵ����
            if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % ����ʽԼ��
            if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % ��ʽԼ��
            obj = obj@func_class(specific);
        end
    end
end
%**************************************************************************
% Reference :   
%   Liang J J, Qu B Y, Suganthan P N. Problem definitions and 
%   evaluation criteria for the CEC 2014 special session and competition 
%   on singl objective real-parameter numerical optimization[J]. 2013.
%**************************************************************************