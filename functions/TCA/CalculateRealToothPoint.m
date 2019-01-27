function RealTP=CalculateRealToothPoint(BasicParameter,ProfileError,mc)
% CalculateRealToothPoint.m             计算带有制造误差的齿面点
% BasicParameter                        摆线轮的基本参数
% ProfileError                          齿廓误差
% mc                                    失配修形系数
% RealTP                                输出实际齿面点

model=mc(3);

zc=BasicParameter(1);
zp=BasicParameter(2);
rp=BasicParameter(3);
a=BasicParameter(4);
rrp=BasicParameter(5);
d_rrp=BasicParameter(6);
d_rp=BasicParameter(7);
deta=BasicParameter(8);
N_PE=length(ProfileError);
rb=a*zp;                          %针轮节圆

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_h=zp/zc;
k1=a*zp/(rp-d_rp);
K=[0,0,1];          %z轴方向单位向量
R_ca=rp+a-rrp;       %摆线轮的齿顶圆半径
Ria=rp-a-rrp;       %摆线轮的齿根圆半径

%       摆线轮方程
syms alfa real
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=mc(4);                                     %-0.002;-0.003;
switch model
    case 2
        %2阶抛物线修形系数及修形参考点
        a2=0.000012345;                      %抛物线修形系数0.000006155669      0.000014815    0.0000098765
        alfa0= 0.056921539864634;                                                          %求得失配修形参考点
        % 2阶抛物线修形摆线轮方程
        equa_c1=[-(rp*sin(alfa+pi)-a*sin(zp*(alfa+pi))+(rrp+b+a2*(sqrt(rp^2+rb^2-2*rp*rb*cos(zc*(alfa+pi)))-sqrt(rp^2+rb^2-2*rp*rb*cos(zc*(alfa0+pi))))^2)*(k1*sin(zp*(alfa+pi))-sin(alfa+pi))/sqrt(1+k1^2-2*k1*cos(zc*(alfa+pi))));
                -(rp*cos(alfa+pi)-a*cos(zp*(alfa+pi))-(rrp+b+a2*(sqrt(rp^2+rb^2-2*rp*rb*cos(zc*(alfa+pi)))-sqrt(rp^2+rb^2-2*rp*rb*cos(zc*(alfa0+pi))))^2)*((-k1)*cos(zp*(alfa+pi))+cos(alfa+pi))/sqrt(1+k1^2-2*k1*cos(zc*(alfa+pi))));
                0;
                1];             %
   case 4   
        %4阶抛物线修形系数及修形参考点
        a2=mc(1);                      %抛物线修形系数0.000006155669     0.00000000182899; 
        alfa0=mc(2);                   %求得失配修形参考点  0.056921539864634; 
        % 4阶抛物线修形摆线轮方程   
        equa_c1=[-(rp*sin(alfa+pi)-a*sin(zp*(alfa+pi))+(rrp+b+a2*(sqrt(rp^2+rb^2-2*rp*rb*cos(zc*(alfa+pi)))-sqrt(rp^2+rb^2-2*rp*rb*cos(zc*(alfa0+pi))))^4)*(k1*sin(zp*(alfa+pi))-sin(alfa+pi))/sqrt(1+k1^2-2*k1*cos(zc*(alfa+pi))));
                -(rp*cos(alfa+pi)-a*cos(zp*(alfa+pi))-(rrp+b+a2*(sqrt(rp^2+rb^2-2*rp*rb*cos(zc*(alfa+pi)))-sqrt(rp^2+rb^2-2*rp*rb*cos(zc*(alfa0+pi))))^4)*((-k1)*cos(zp*(alfa+pi))+cos(alfa+pi))/sqrt(1+k1^2-2*k1*cos(zc*(alfa+pi))));
                0;
                1]; 
end
% %齿根在上
% equa_c1=[-(((rp-d_rp)-(rrp+d_rrp)/sqrt(1+k1^2-2*k1*cos(zc*alfa)))*sin((1-i_h)*zc*alfa-deta)+a/(rp-d_rp)*(rp-d_rp-zp*(rrp+d_rrp)/sqrt(1+k1^2-2*k1*cos(zc*alfa)))*sin(i_h*zc*alfa+deta));
%          (((rp-d_rp)-(rrp+d_rrp)/sqrt(1+k1^2-2*k1*cos(zc*alfa)))*cos((1-i_h)*zc*alfa-deta)-a/(rp-d_rp)*(rp-d_rp-zp*(rrp+d_rrp)/sqrt(1+k1^2-2*k1*cos(zc*alfa)))*cos(i_h*zc*alfa+deta));
%           0;
%           1]; 
%对摆线轮1齿廓求单位法向量
d_equa_c1_alfa=diff(equa_c1,'alfa',1);           %切矢
d1=d_equa_c1_alfa(1:3);
N_c1=cross(d1,K);                                %叉乘，得到法向量
M0=sqrt(N_c1(1)^2+N_c1(2)^2+N_c1(3)^2);          %法向量的模
n_c1=N_c1/M0;                                    %摆线轮上接触点处的单位法向量
alpha=linspace(0,pi/zc,N_PE);
x_sheji=zeros(N_PE,1);% 计算设计齿面点
y_sheji=zeros(N_PE,1);
for i=1:N_PE
    alfa0=alpha(i);
    xs0=subs(equa_c1(1),'alfa',alfa0);
    ys0=subs(equa_c1(2),'alfa',alfa0);
    x_sheji(i)=xs0;
    y_sheji(i)=ys0;
end
% plot(x_sheji,y_sheji,'b')
% hold on

x_real=zeros(N_PE,1);
y_real=zeros(N_PE,1);
for i=1:N_PE
    x0_real=equa_c1(1)-n_c1(1)*ProfileError(i);
    y0_real=equa_c1(2)-n_c1(2)*ProfileError(i);
    alfa0=alpha(i);
   x0=subs(x0_real,'alfa',alfa0);
   y0=subs(y0_real,'alfa',alfa0);
   x_real(i)=x0;
   y_real(i)=y0;
end
% plot(x_real,y_real,'r')
% hold off
RealTP=[x_real,y_real];

