function [x,y]=equation_fit(i,p,U,DPoint)
% BasisFuns.m       计算出基函数的值N_i-k(u)、N_i-k+1(u)、...、N_i(u)
% i                 节点区间的下标
% u                 给定的u值
% p                 拟合曲线的次数
% U                 节点矢量
% CPoint            控制顶点

syms u
right=sym(zeros(1,p));
left=sym(zeros(1,p));
N=sym(zeros(1,p+1));
N(1)=1;
for j=2:p+1
    left(j)=u-U(i+1+1+1-j);
    right(j)=U(i+j)-u;
    saved=0;
    for r=1:j-1
        temp=N(r)/(right(r+1)+left(j-r+1));
        N(r)=saved+right(r+1)*temp;
        saved=left(j-r+1)*temp;
    end
    N(j)=saved;
end
I=i;
Curve0(1,:)=N(1).*DPoint(I+1-3,:)+N(2).*DPoint(I+1-2,:)+N(3).*DPoint(I-1+1,:)+N(4).*DPoint(I+1,:);
x=Curve0(1);
y=Curve0(2);
