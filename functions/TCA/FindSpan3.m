function N=FindSpan3(ObjParam,ScopeParam)
% FindSpan3.m           在一组数组中找到目标值所在的区间,与findspan2的区别，有待修改
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

low=1;high=N_ScopeParam;
mid=fix((low+high)/2);
% if ObjParam <= ScopeParam(1)&&ObjParam >= ScopeParam(2)
%     N=2;
%     return;
% end
% if ObjParam >= ScopeParam(end)&&ObjParam <= ScopeParam(end-1)
%     N=N_ScopeParam;
%     return;
% end

while ObjParam < ScopeParam(mid)||ObjParam >= ScopeParam(mid+1)
   if ObjParam < ScopeParam(mid)
        high=mid;
   else
        low=mid;
   end
   mid=fix((low+high)/2);
end
N=mid+1;