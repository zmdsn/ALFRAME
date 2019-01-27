%利用牛顿迭代法求解方程组的解

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%给出实际齿面点，利用b样条拟合方法拟合曲线
RealPoint=load('realtoothpoint2.txt');
[ControlPoint,NodeVector]=BYT_Fitting(RealPoint,3);
N_DPoint=length(ControlPoint(:,1));
U=NodeVector;                                    %将前面计算的节点矢量赋给U
DPoint=ControlPoint;                             %将前面计算的控制顶点矢量赋给DPoint

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%先建立求解误差的方程
%理论摆线轮部分
%%%%%%%%%   摆线轮主要设计参数    %%%%%%%%%
zc=41;                                          %摆线轮齿数%
zp=42;                                          %针轮齿数%
rp=70;                                          %针齿分布圆半径%
a=1;                                            %偏心距%
rrp=2.5;                                        %针齿半径
d_rrp=0;                                        %等距修形量(正+表示砂轮半径增加)
d_rp=0;                                         %移距修形量（+表示针齿分布半径减小）
deta=0;                                         %转角修形量
%
i_h=zp/zc;
k10=a*zp/rp;                                    %理论短幅系数
k1=a*zp/(rp-d_rp);                              %修形之后短幅系数
R_ca=(rp-d_rp)+a-(rrp+d_rrp);                   %摆线轮的齿顶圆半径
R_ia=(rp-d_rp)-a-(rrp+d_rrp);                   %摆线轮的齿根圆半径
K=[0,0,1];                                      %z轴方向单位向量

syms alfa x_fit y_fit miu real
% alfa   为摆线轮c1的齿廓方程参数
% x_fit  为拟合方程的横坐标 
% y_fit  为拟合方程的纵坐标
% miu    为齿面误差量

%%%%%%%      摆线轮方程    %%%%%%%%%
equa_c1=[-(((rp-d_rp)-(rrp+d_rrp)/sqrt(1+k1^2-2*k1*cos(zc*alfa)))*sin((1-i_h)*zc*alfa-deta)+a/(rp-d_rp)*(rp-d_rp-zp*(rrp+d_rrp)/sqrt(1+k1^2-2*k1*cos(zc*alfa)))*sin(i_h*zc*alfa+deta));
         (((rp-d_rp)-(rrp+d_rrp)/sqrt(1+k1^2-2*k1*cos(zc*alfa)))*cos((1-i_h)*zc*alfa-deta)-a/(rp-d_rp)*(rp-d_rp-zp*(rrp+d_rrp)/sqrt(1+k1^2-2*k1*cos(zc*alfa)))*cos(i_h*zc*alfa+deta));
          0;
          1]; 
%对摆线轮1齿廓求单位法向量
d_equa_c1_alfa=diff(equa_c1,'alfa',1);           %切矢
d1=d_equa_c1_alfa(1:3);
N_c1=cross(d1,K);                                %叉乘，得到法向量
M0=sqrt(N_c1(1)^2+N_c1(2)^2+N_c1(3)^2);          %法向量的模
n_c1=N_c1/M0;                                    %摆线轮上接触点处的单位法向量
x0=equa_c1(1);
y0=equa_c1(2);

F1=x0+miu*n_c1(1)-x_fit;
F2=y0+miu*n_c1(2)-y_fit;
dj=10^(-6);

miu0=0.0241;
u0=0.4773;
C=CurvePoint(u0,3,U,DPoint);
f1_k=subs(subs(subs(F1,'alfa',0.02766985),'miu',miu0),'x_fit',C(1));
f2_k=subs(subs(subs(F2,'alfa',0.02766985),'miu',miu0),'y_fit',C(2));
      
 


Jcb()

