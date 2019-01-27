function [ f_rp1,f_np1 ] = Surfacepoint_Pin(BasicParameter1,beta0,fai1)
% Surfacepoint_Pin.m        针齿再固定坐标系中的齿面点
% BasicParameter1           输入摆线针轮基本参数
% beta0                     输入针齿参数
% fai1                      输入针轮转角
% f_rp1                     输出针轮齿面点
% f_np1                     输出针齿法向单位矢量

K=[0,0,1];
rp=BasicParameter1(3);
rrp=BasicParameter1(5);


% 针轮到固定坐标系的坐标转换%
M_fp1=[cos(fai1),-sin(fai1),0,0;
       sin(fai1),cos(fai1),0,0;
       0,0,1,0;
       0,0,0,1];
syms beta1 real
% 针齿方程
equa_r1=[-rrp*sin(beta1);rp-rrp*cos(beta1);0;1];

%对针齿方程求单位法向量
d_equa_r1_beta1=diff(equa_r1,'beta1',1);                %针齿求切矢
d2=d_equa_r1_beta1(1:3)';
N_p1=cross(K,d2);                                       %N_p1的值为复数，可能有问题，后面可能需要调整
M1=sqrt(N_p1(1)^2+N_p1(2)^2+N_p1(3)^2);                 %法向量的模
n_p1=N_p1/M1;                                           %针轮上接触点处的单位法向量
% rx1(1)=subs(equa_r1(1),'beta1',beta0);                         %求位置矢量的数值解
% ry1(2)=subs(equa_r1(2),'beta1',beta0);                         %求位置矢量的数值解
% n1(1)=subs(n_p1(1),'beta1',beta0);                            %求法向矢量的数值解
% n1(2)=subs(n_p1(2),'beta1',beta0);                            %求法向矢量的数值解
r1=double(subs(equa_r1,'beta1',beta0));                         %求位置矢量的数值解
n1=double(subs(n_p1,'beta1',beta0));                            %求法向矢量的数值解


f_rp1=M_fp1*r1;
f_np1=M_fp1(1:3,1:3)*n1';

end

