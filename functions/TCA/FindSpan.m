function I=FindSpan(N,p,u,U)
% FindSpan.m            用二分法确定节点下标
% N                     节点矢量个数
% p                     拟合曲线次数
% u                     给定的节点值
% U                     节点矢量
% I                     输出的节点下标值

n_DP=N-p-2;                %控制顶点的最大下标,n=m-p-1,m为节点最大下标或者n=numel(U)-p-2
if u==U(n_DP+2)|| abs(u-U(n_DP+2))<10^(-6)          %特殊情况：u如果是定义域内的最后一个值，直接返回n
    I=n_DP;
    return
end
low=p+1;high=n_DP+2;
mid=fix((low+high)/2);
while u<U(mid)||u>=U(mid+1)
   if u<U(mid)
        high=mid;
   else
        low=mid;
   end
   mid=fix((low+high)/2);
end
I=mid-1;                                            %返回的结果和书本上面的结果相同



    
