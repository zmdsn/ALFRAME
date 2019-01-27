classdef cec14_f22  < func_class 
    methods       
        %% �����㺯��
        function  [f, g, h] = object_function(obj, x , specific )
            % ע������xΪ�б���
            g = zeros(obj.Fneq,1);
            h = zeros(obj.Feq,1);
            f = cec14_func(x,22);
        end
        %% ���캯��
        function obj = cec14_f22(specific)
            % �������½��Լ�ά��
            % N = 5
            % p = [0.1, 0.2, 0.2, 0.2, 0.3]
            % g1 : Katsuura Function f10
            % g2 : HappyCat Function f11
            % g3: Expanded Griewank��s plus Rosenbrock��s Function f13
            % g4: Modified Schwefel��s Function f9
            % g5: Ackley��s Function f5
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