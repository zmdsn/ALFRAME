function [ a ] = randnSm( n,m ,k)
% rand select m from 1:n
    if nargin < 3
        k = n;
    end
    if n < m
       error('please check your  input : m>n') 
    end
    a = randperm(n);
    
    if find(a==k) <= m
        a = setdiff( a(1:m+1) , k);
    else
        a = a(1:m);
    end
end

