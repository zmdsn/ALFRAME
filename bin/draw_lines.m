function draw_lines(fun_data , nlines, Run , func_num ,lengs_str, type ,folds)
% 绘图函数,用于绘制一条线
    HowMuchStar = 20;  % how large the interval is to plot one dot
    figure;
    set(0,'defaultfigurecolor','w') 
    if nargin < 6
        type = 3;
    end
    % type 取值1,2,3分别代表不同的做表取值
    for j = 1:nlines
      draw_one_line( j ,fun_data{j} ,HowMuchStar ,type, Run) 
    end
    legend(lengs_str);
    xlabel('\fontsize{15}Function Evaluations');
    ylabel('\fontsize{15}Average Function Values');
    title(['\fontsize{14}Comparison among algorithms on f',num2str(func_num), ' with ', num2str(Run), ' runs'])
    if ~exist('folds','var') 
        if ~exist('./figure/') 
            mkdir('./figure/');         % 若不存在，在当前目录中产生一个子目录‘Figure’
        end
        figure_name = ['./figure/fun'];
    else
        folds_name = [folds 'figure/'];
        if ~exist(folds_name,'dir') 
            mkdir(folds_name);         % 若不存在，在当前目录中产生一个子目录‘Figure’
        end
        figure_name = [folds_name 'fun'];
    end
    saveas(gcf,[figure_name,num2str(func_num)],'fig');
end

function draw_one_line(order , line_data ,HowMuchStar, type, Run)  
    A = line_data(:,2)./Run;
    line_color = {'r','m','b','k','g','y','c','r','m','b','k','g','y','c','r','m','b','k','g','y','c'};
    line_style = {'*','o','>','+','diamond','x' ,'<', 'square','v','pentagram','hexagram','*','o','>','+','*','o','>','+','diamond','x' ,'<', 'square','v'};
    switch type
        case 1
          p1 =semilogy(line_data(:,1), A,line_color{order});hold on;
          p11=semilogy(line_data(1:floor(length(A)/HowMuchStar):end,1), A(1:floor(length(A)/HowMuchStar):end),[line_color{order} line_style{order}],  'MarkerSize',8);hold on;
          p12=semilogy(line_data(1,1), A(1),line_color{order},'Marker',line_style{order},'MarkerSize',8);
          set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
          set(get(get(p11,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
          set(p1,'LineWidth',2.0);
          hold on
        case 2 
          p1 =semilogx(line_data(:,1), A,line_color{order});hold on;
          p11=semilogx(line_data(1:floor(length(A)/HowMuchStar):end,1), A(1:floor(length(A)/HowMuchStar):end),[line_color{order} line_style{order}],  'MarkerSize',8);hold on;
          p12=semilogx(line_data(1,1), A(1),line_color{order},'Marker',line_style{order}, 'MarkerSize',8);
          set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
          set(get(get(p11,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
          set(p1,'LineWidth',2.0);
          hold on
        otherwise
          % 显示曲线    
          p1 = plot(line_data(:,1), A,line_color{order});hold on;
          % 显示曲线上的点
%           p11=plot(line_data(1:floor(length(A)/HowMuchStar):end,1), A(1:floor(length(A)/HowMuchStar):end),[line_color{order} line_style{order}], 'MarkerSize',8);hold on;
          p11=plot(line_data(1:floor(length(A)/HowMuchStar):end,1), A(1:floor(length(A)/HowMuchStar):end),[line_color{order} line_style{order}], 'MarkerFaceColor',line_color{order}, 'MarkerSize',8);hold on;
          % 显示图例点
          p12=plot(line_data(1,1), A(1),line_color{order},'Marker',line_style{order}, 'MarkerFaceColor',line_color{order}, 'MarkerSize',8);
%           p12=plot(line_data(1,1), A(1),line_color{order},'Marker',line_style{order},  'MarkerSize',8);
          % 隐藏p1,p11的图例
          set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
          set(get(get(p11,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
          set(p1,'LineWidth',2.0);
          hold on
    end
end