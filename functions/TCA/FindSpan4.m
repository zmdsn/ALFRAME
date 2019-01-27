function N=FindSpan4(ObjParam,ScopeParam)
% FindSpan3.m           在一组数组中找到目标值所在的区间,与findspan2的区别有待研究
% ObjParam              输入目标值
% ScopeParam            输入范围数组
% N                     输出比目标值大，但最接近的一个数的序号


N_ScopeParam=length(ScopeParam);
% 首先判断是否超出搜索范围
Max_sp = max(ScopeParam);
Min_sp = min(ScopeParam);
if ObjParam >= Max_sp||ObjParam < Min_sp 
      N=2;
      return;
end
%如果不超出搜索范围继续求解
%对两端进行判断
if ObjParam <= ScopeParam(1)&&ObjParam >= ScopeParam(2)
    N=2;
    return;
end
if ObjParam >= ScopeParam(end)&&ObjParam <= ScopeParam(end-1)
    N=N_ScopeParam;
    return;
end

% 对中间值进行判断
high=1;low=N_ScopeParam;
mid=fix((low+high)/2);

while ObjParam > ScopeParam(mid)||ObjParam <= ScopeParam(mid+1)
   if ObjParam > ScopeParam(mid)
        low=mid;
   else
        high=mid;
   end
   mid=fix((low+high)/2);
end
N=mid+1;