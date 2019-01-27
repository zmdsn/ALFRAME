function N=FindSpan2(ObjParam,ScopeParam)
% FindSpan2.m           在一组数组中找到目标值所在的区间
% ObjParam              输入目标值
% ScopeParam            输入范围数组
% N                     输出比目标值大，但最接近的一个数的序号

N_ScopeParam=length(ScopeParam);
low=1;high=N_ScopeParam;
mid=fix((low+high)/2);
while ObjParam>ScopeParam(mid)||ObjParam<=ScopeParam(mid+1)
   if ObjParam>ScopeParam(mid)
        high=mid;
   else
        low=mid;
   end
   mid=fix((low+high)/2);
end
N=mid;