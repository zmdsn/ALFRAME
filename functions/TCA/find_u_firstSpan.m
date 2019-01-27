function U_First_No=find_u_firstSpan(R_half,U,DPoint,p)
% find_u_first              计算出分度圆半径处的u值
% R_half                    输入分度圆半径
% U                         输入节点矢量
% DPoint                    输入控制顶点
% p                         拟合样条次数
% U_First_No                输出初始参考点的u值坐在区间的节点下标

N_DPoint=length(DPoint);
SR_half=R_half^2;
m=1;
Pt=zeros(1,2);
for i=p+1:N_DPoint+1
    I=i;
    u1=U(I);
%     u2=U(I+1);
    P_fit1=Point_of_fitting(u1,U,DPoint);
    Pt(m,:)=P_fit1;
%     P_fit2=Point_of_fitting(u2,U,DPoint);
%     SP1=P_fit1(1)^2+P_fit1(2)^2;
%     SP2=P_fit2(1)^2+P_fit1(2)^2;
    
%     if SR_half>SP2 && SR_half<= SP1
%        break 
%     end
    m=m+1;
end
SPt=Pt(:,1).^2+Pt(:,2).^2;
No_SPt=FindSpan2(SR_half,SPt);
I=No_SPt+3-1;   %输出初始参考点的u值坐在区间下标          
U_First_No=I;

