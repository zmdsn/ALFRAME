function N=BasisFuns(i,u,p,U)
% BasisFuns.m       计算出基函数的值N_i-k(u)、N_i-k+1(u)、...、N_i(u)
% i                 节点区间的下标
% u                 给定的u值
% p                 拟合曲线的次数
% U                 节点矢量

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