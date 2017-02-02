function [R] = start_product(V,x,b,k)
    n = size(V,1);
    R = zeros(n,b);
    for i=1:k
        R = R + x(i)*V(:,(i-1)*b+1:i*b);
    end
    
    
    