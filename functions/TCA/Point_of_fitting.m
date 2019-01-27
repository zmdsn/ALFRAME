function P_fit=Point_of_fitting(u,U,DPoint)
% Point_of_fitting              计算拟合齿面点
% u                             输入拟合参数u
% U                             输入节点矢量
% DPoint                        输入控制顶点
% P_fit                         输出计算得到的拟合齿面点

N_U=length(U);
I=FindSpan(N_U,3,u,U);            %第I个节点在matlab中是第I+1个，所以在计算曲线时，都在I上加了1
N=BasisFuns(I,u,3,U);
C=N(1).*DPoint(I+1-3,:)+N(2).*DPoint(I+1-2,:)+N(3).*DPoint(I-1+1,:)+N(4).*DPoint(I+1,:);
P_fit=C;