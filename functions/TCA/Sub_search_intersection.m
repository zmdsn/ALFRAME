function [X1,f1]=Sub_search_intersection(A1,A2)
% Sub_search_intersection.m     对离散的数据点拟合，寻找交点
% A1                            为输入曲线1参数
% A2                            为输入曲线2参数
% X1                            为输出交点处的横坐标值
% f1                            为输出交点处的纵坐标值

N_A1=length(A1(:,1));
for i=1:N_A1
N=FindSpan4(A1(i,1),A2(:,1));        %二分法求出给定点所在范围
a1=A1(i,2)-A2(N-1,2);
a2=A1(i,2)-A2(N,2);
if abs(a1)<0.5||abs(a2)<0.5
    N1=i;
    N2=N;
    break;
end
end

p1 = [A1(N1-2,:);A1(N1-1,:);A1(N1,:);A1(N1+1,:)];            % 以两个最近点的位置，分别取三个点，用来3次多项式拟合
p2 = [A2(N2-2,:);A2(N2-1,:);A2(N2,:);A2(N2+1,:)];
  
x1= p1(:,1)';
y1=p1(:,2)';
x2= p2(:,1)';
y2=p2(:,2)';
  

  P1=polyfit(x1,y1,3);                              % 第一条曲线的3次多项式拟合
  P2=polyfit(x2,y2,3);                              % 第二条曲线的3次多项式拟合
  
%   绘制拟合曲线
  xx1=x1(1):0.01:x1(end);
  yy1=polyval(P1,xx1);
  xx2=x2(1):0.01:x2(end);
  yy2=polyval(P2,xx2); 
  
%两曲线做差，求两曲线交点
P12=P1-P2 ;
xx3=x1(1):0.01:x2(end);
% xx3=-5:0.5:150;
yy3=polyval(P12,xx3);


%   绘制图像
%绘制原数据图像
  plot(x1,y1)
  hold on 
  plot(x2,y2)
%绘制拟合数据图像  
  plot(xx1,yy1,'b')
%   hold on 
  plot(xx2,yy2,'r')
  plot(xx3,yy3,'k')
  
  hold off
  grid on 
  
% P12的多项式方程表达式
syms x
F12=P12(1)*x^3+P12(2)*x^2+P12(3)*x+P12(4);
  
%建立方程函数文件
fid=fopen('equation_of_search_point.m','wt+');
fprintf(fid,'%s\n','function f=equation_of_search_point(x)');    
fprintf(fid,'%s','f=');
fprintf(fid,'%s\n',strcat(char(F12),';'));
fclose(fid);
% 求解多项式方程
X0=x1(2);
[x,fval]=fzero('equation_of_search_point',X0);
%   r=roots(P12)
X1=x;                                               %输出两曲线交点处的横坐标
f1=polyval(P1,X1);                                  %输出两曲线交点处的横坐标
      






