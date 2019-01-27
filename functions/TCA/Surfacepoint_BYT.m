function [ f_rc1,f_nc1 ] = Surfacepoint_BYT( u0,p,U,DPoint,U_d,Q1,fai2,a)
%Surfacepoint_BYT.m        求解出拟合齿面在转化轮系中的点 
% u0                       拟合参数u
% p                        拟合次数
% U                        节点矢量
% DPoint                   控制顶点
% U_d                      法向节点矢量
% Q1                       法向控制顶点
% fai2                     摆线轮转过角度
% a                        偏心距
% f_rc1                    固定坐标系中摆线轮坐标
% f_nc1                    固定坐标系中摆线轮法向单位矢量

% 坐标转换矩阵
M_fc1=[cos(fai2),-sin(fai2),0,0;
    sin(fai2),cos(fai2),0,-a;
    0,0,1,0;
    0,0,0,1];

[C,n_DC]=CurvePoint(u0,p,U,DPoint,U_d,Q1);
rc_point=[C,0,1]';%摆线轮齐次坐标


f_rc1=M_fc1*rc_point;
f_nc1=M_fc1(1:3,1:3)*n_DC';
end

