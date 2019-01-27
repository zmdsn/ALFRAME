classdef func_class 
% ����:
    % Xmin Ĭ���Ա����½�
    % Xmin Ĭ���Ա����Ͻ�
    % Xdim �Ա���ά��
    % Fmin Ŀ�꺯����Сֵ
    % Tex  Ŀ�꺯��Latex����
% ����:
    % func_class(specific)                          ���캯��
    % object_function(obj , x ,specific)               Ŀ�꺯��
    % fun_figure                                    1 1ά����,2 ��ά����͵ȸ���?3 �ȸ���
    % initialize(obj, sw_size )                     ��ʼ����Ⱥ;
    % draw_running(obj,swarm)                       �㷨���в��պ���
    % draw_road(obj,citys,road,specific)            ����·��ͼ(TSP)
    % multi_draw(obj,f,M)                           ���ƶ�Ŀ�꺯��,Ŀ��������?    % is_one2n(obj,x)                               �ж��Ƿ�Ϊ1:n�����е�����
    %% �������?    
    properties
        Xmin = [];  % �½�
        Xmax = [];  % �Ͻ�
        Xdim = 30;  % ����ά��
        Tex = '';   % ����latex����?        Fmin = 0;   % ������Сֵ
        Nobj = 1;   % Ŀ�꺯��ֵ����
        Fneq = 0;   % ����ʽԼ��
        Feq  = 0;   % ��ʽԼ��
        Fmin = 0;
    end %properties
    methods(Abstract)
       %% �����㺯��
        [f, g, h] = object_function(obj,x,specific);
    end
    methods
       %% ���캯��
        function obj = func_class(specific)
            % initialize variable ; Xdim
            if isfield(specific,'Xdim') 
                if length(specific.Xdim) == 1
                    obj.Xdim = specific.Xdim;
                else
                    error('Please enter the correct a parameter : Xdim.Please be sure the dimension of Xdim is 1');
                end
            else
                warning(['The dimension of  variables is set as the default Xdim= 30.'...
                       'If it is not your original intention, please check the parameter : specific.Xdim']);
                obj.Xdim = 30;
            end
            % initialize variable ; Xmin
            if isfield(specific,'Xmin')
                if length(specific.Xmin) == 1
                    obj.Xmin = specific.Xmin*ones(obj.Xdim,1);
                else
                    obj.Xmin = specific.Xmin;
                end
            else
                warning(['The dimension of  variables is set as the default Xmin= -100.'...
                       'If it is not your original intention, please check the parameter : specific.Xmin']);
                obj.Xmin = -100*ones(obj.Xdim,1);
            end
            % initialize variable ; Xmax
            if isfield(specific,'Xmax')
                if length(specific.Xmax) == 1
                    obj.Xmax = specific.Xmax*ones(obj.Xdim,1);
                else
                    obj.Xmax = specific.Xmax;
                end
            else
                warning(['The dimension of  variables is set as the default Xmax= -100.'...
                       'If it is not your original intention, please check the parameter : specific.Xmax']);
                obj.Xmax = -100*ones(obj.Xdim,1);
            end
            % initialize variable ; Fmin
            if isfield(specific,'Fmin') 
                if length(specific.Fmin) == 1
                    obj.Fmin = specific.Fmin;
                else
                    error('Please enter the correct parameters : Fmin');
                end
            else
                warning(['The dimension of  variables is set as the default Fmin= -100.'...
                       'If it is not your original intention, please check the parameter : specific.Fmin']);
                obj.Fmin = 0;
            end
            % initialize variable ; Tex
            if isfield(specific,'Tex') 
                obj.Tex = specific.Tex;
            else
                obj.Tex = 'lazy man without input any Latex code';
            end
            % initialize variable ; Nobj
            if isfield(specific,'Nobj') 
                obj.Nobj = specific.Nobj;
            else
                obj.Nobj = 1;
            end
            % initialize variable ; Fneq
            if isfield(specific,'Fneq') 
                obj.Fneq = specific.Fneq;
            else
                obj.Fneq = 0;
            end
            % initialize variable ; Feq
            if isfield(specific,'Feq') 
                obj.Feq = specific.Feq;
            else
                obj.Feq = 0;
            end
        end % end func_class
        
        % calculate the value of the objective function
        function [f, g, h] = obj_value(obj, x , specific )
            % ע������xΪ�б���
            % specific Ԥ��ֵ,����ʹ�ýṹ����,���ڴ��ݶ���Ĳ���?            if ~exist('specific','var')
                specific = [];  % �����ʼ��?            end
            swarmsize = size( x , 2 );
            if swarmsize > 1
            % ����һ����Ⱥ�ĺ���ֵ
				f = zeros( obj.Nobj , swarmsize ); 
                g = zeros( obj.Fneq , swarmsize ); 
                h = zeros( obj.Feq  , swarmsize ); 
                for i = 1:swarmsize
                    [f(:,i), g(:,i) ,h(:,i)]  = obj.object_function(x(:,i),specific);
                end
            else
            % ����һ�������ĺ���ֵ
                [f, g, h] = obj.object_function(x,specific);
%                 f = obj.object_function(x,specific);
            end
        end % end obj_value
        
       %% ������m���Ա���
        function init_matrices = initialize(obj, sw_size )
            % ������һ����Ⱥ����һ����,Ĭ����Ⱥ�����ÿһ�д��һ����
            switch nargin
                case 1
                    init_matrices = rand(obj.Xdim , 1 ).*(obj.Xmax - obj.Xmin) + obj.Xmin;
                otherwise
                    rep_Xmax = repmat(obj.Xmax,[1 sw_size]);
                    rep_Xmin = repmat(obj.Xmin,[1 sw_size]);
                    init_matrices = rand( obj.Xdim , sw_size).*(rep_Xmax - rep_Xmin) + rep_Xmin;
            end
        end % end initialize
        
       %% ��ͼ��ʾ����
        function  [] = fun_figure(obj,type)
            % ��ͼ����
            h = (obj.Xmax - obj.Xmin)/150;%���?  
            x1 = [obj.Xmin:h(1):obj.Xmax];
            y1 = [obj.Xmin:h(1):obj.Xmax];
            [x,y] = meshgrid(x1,y1);
            switch type
                case 3
                    % �Ա���Ϊ��ά,��ʾΪ��ά����͵ȸ���?                    [x,y] = meshgrid(obj.Xmin(1):h(1):obj.Xmax(1),obj.Xmin(2):h(2):obj.Xmax(2));
                    [a,b] = size(x);
                    z = zeros(size(x));
                    for i=1:a
                        for j=1:b
                            z(i,j) = obj.obj_value([x(i,j) y(i,j)]');
                        end
                    end
                    [c,h]=contour(x,y,z,60);
                case 2
                    % �Ա���Ϊ��ά,��ʾΪ��ά����͵ȸ���?                    [x,y] = meshgrid(obj.Xmin(1):h(1):obj.Xmax(1),obj.Xmin(2):h(2):obj.Xmax(2));
                    [a,b] = size(x);
                    z = zeros(size(x));
                    for i=1:a
                        for j=1:b
                            z(i,j) = obj.obj_value([x(i,j) y(i,j)]');
                        end
                    end
                    meshc(x,y,z);% or surfc(x,y,z);
                case 1
                    x = [obj.Xmin:h(1):obj.Xmax];
                    fun_matrices = obj.obj_value( x );
                    plot(x,fun_matrices);
                otherwise
                    error('�ú����Ա���������ά�޷�չʾ');
            end
        end % end fun_figure
        
        % �����㷨�е�״̬
        function [] = draw_running(obj,swarm,fname)   
            hold off;
            scatter(swarm(1,:),swarm(2,:));
            hold on;
            obj.fun_figure(3);
            pause(.1);
            if exist('fname','var')
                f = getframe(gcf);  
                imind = frame2im(f);
                [im,map] = rgb2ind(imind,256);
                gifname = ['GIFS/' fname '.gif'];
                if ~exist(gifname,'file')
                    imwrite(im,map,gifname,'gif','loopcount',inf)
                else
                    imwrite(im,map,gifname,'gif','writemode','append')
                end
            end
        end
        
        % ����TSP·��
        function [] = draw_road(obj,citys,road,specific)
            plot( citys(road,1),citys(road,2),'o-');
            grid on
            for i = 1:size(citys,1)
                text(citys(i,1),citys(i,2),['   ' num2str(i)]);
            end
            text(citys(road(1),1),citys(road(1),2),'        start','Color','red');
            text(citys(road(end-1),1),citys(road(end-1),2),'        end','Color','red');
            
            if isfield(specific,'xlabel')
                xlabel(specific.xlabel);
            else
                xlabel('�����');
            end
            if isfield(specific,'ylabel')
                ylabel(specific.ylabel);
            else
                ylabel('�����');
            end
            if isfield(specific,'title')
                title(specific.title);
            else
                title('This is title');
            end
        end % end draw_road
        
        % ���ƶ�Ŀ��ĺ�����
        function [] = multi_draw(obj,f,M)
            % f Ϊ�����������Ŀ�꺯������?�������ĸ���
            figure;
            switch M
                case 1
                    warning('This function is used for testing multiobjective functions.Please check again ');
                case 2
                    plot(f(1,:),f(2,:),'*');
                case 3
                    plot3(f(1,:),f(2,:),f(3,:),'*');  
                otherwise                    
%                     set(fid,'paperposition',[0.3,6,20,15])
%                     get(fid,'paperposition')
                    set(0,'defaultfigurecolor','w')   
%                     set(gca,'position',[0.01,0.01,0.9,0.9]);
                    H = zeros(M);
                    for i = 1:M
                        for j = 1:M
                            H(M*(i-1)+j) = subplot(M,M,M*(i-1)+j);
                            if i == j
                                text(0.45,0.5,['f' num2str(i)]);
                                set(gca,'xtick',[],'xticklabel',[])
                                set(gca,'ytick',[],'yticklabel',[])
                            else
                                plot(f(i,:),f(j,:),'*');
                                set(gca,'xtick',[],'xticklabel',[])
                                set(gca,'ytick',[],'yticklabel',[])
                            end
                            PPP=get(H(M*(i-1)+j),'pos');%��NN����ͼ�ĵ�ǰλ��
                            PPP(1)=PPP(1)-0.03;%���ұ���չ0.04
                            PPP(2)=PPP(2)-0.03;%���ұ���չ0.04
                            PPP(3)=PPP(3)+0.02;%���ұ���չ0.04
                            PPP(4)=PPP(4)+0.03;%���Ϸ���չ0.03
                            set(H(M*(i-1)+j),'pos',PPP)%����µı߽����á�?                        
                        end
                    end
            end
        end % end multi_draw
                
        % �ж�һ������Ƿ��?��n����������
        function result = is_one2n(obj, x)
            n = size(x,1) * size(x,2);
            if size(x,1) ~= 1 && size(x,2) ~= 1
               error('dimension error');
            end
            x = reshape(x,n,1);
            array = [1:n]';
            result = sum((sort(x) - array) == zeros(n,1)) == n;
        end
        
    end %methods
end