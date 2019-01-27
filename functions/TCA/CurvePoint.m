function [C,n_DC]=CurvePoint(u,p,U,DPoint,U_d,Q1)
% CurvePoint.m              给出曲线参数u值，计算出拟合曲线上的点
% n                         是控制顶点的最大下标或者是控制顶点个数-1
% u                         拟合曲线上的参数值u
% DPoint                    控制顶点
% U                         节点矢量
% C                         输出结果 
% n_DC                      拟合曲线的法向矢量
N_U=length(U);
K=[0,0,1];
I=FindSpan(N_U,p,u,U);            %第I个节点在matlab中是第I+1个，所以在计算曲线时，都在I上加了1
N=BasisFuns(I,u,p,U);
C=N(1).*DPoint(I+1-3,:)+N(2).*DPoint(I+1-2,:)+N(3).*DPoint(I-1+1,:)+N(4).*DPoint(I+1,:);
N_d=BasisFuns(I-1,u,2,U_d);%因为导函数是比原拟合曲线低一阶的曲线拟合，所以节点矢量下标要减1
DC=N_d(1).*Q1((I-1)-2+1,:)+N_d(2).*Q1((I-1)-1+1,:)+N_d(3).*Q1(I-1+1,:);
DC=[DC,0];
N_DC=cross(DC,K);
n_DC=N_DC/sqrt((N_DC(1)^2+N_DC(2)^2));






