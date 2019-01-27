function [ControlPoint,NodeVector]=BYT_Fitting(InPar,k)
% BYT_Fitting.m       该函数利用B样条拟合方法求出控制顶点和节点矢量
% InPar               输入需要拟合的数据点（即实际齿面点）
% k                   需要拟合次数
% ControlPoint        输出控制顶点
% NodeVector          输出节点矢量

[n,m]=size(InPar);          %n是数据点数目，m为数据点维数
%-----------弦长参数化--------------------------------------
L=0;
for i=1:n-1
    L=L+sqrt((InPar(i+1,1)-InPar(i,1))^2+(InPar(i+1,2)-InPar(i,2))^2);
end;
u_k(n)=0;
u_k(1)=0;
for i=1:n-1
    u_k(i+1)=u_k(i)+sqrt((InPar(i+1,1)-InPar(i,1))^2+(InPar(i+1,2)-InPar(i,2))^2)/L;
end;
%用平均值法求节点矢量
U(n+k+1)=0;                                 %此处和书本上（n+k+2）不同，是因为matlab的下标从1开始的
for i=1:(n-1)-3
 U(i+k+1)=1/k*(u_k(i+1)+u_k(i+2)+u_k(i+3)); %此处+1是因为书本上的节点是从零开始，+1方便对应
end
U(n+4-3:n+4)=[1,1,1,1];
%至此参数化已经全部处理完毕
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%构造基函数矩阵
Nb=zeros(n,n);
Nb(1,1)=1;
Nb(n,n)=1;
N_U=length(U);
for i=2:n-1
    uu=u_k(i);
    I =FindSpan(N_U,3,uu,U);
    N =BasisFuns(I,uu,3,U);
    I_1=I-3+1;
    Nb(i,I_1:I_1+3)=N;
end
Va=InPar;
DPoint=Nb\Va;
% plot(InPar(:,1),InPar(:,2),'ro')
% hold on
% plot(DPoint(:,1),DPoint(:,2),'g*')
ControlPoint=DPoint;
NodeVector=U;



