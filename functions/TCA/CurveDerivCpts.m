function Q1=CurveDerivCpts(n,k,U,P)
% CurveDerivCpts.m      计算导曲线的控制顶点
% n                     原来控制顶点的最大下标
% k                     拟合曲线次数
% U                     曲线拟合的节点矢量
% P                     控制顶点
% Q1                    求得导数的控制顶点

Q1=zeros(n,2);
tmp=k;
for i=1:n
    Q1(i,:)=tmp.*((P(i+1,:)-P(i,:))./(U(i+k+1)-U(i+1)));
end

