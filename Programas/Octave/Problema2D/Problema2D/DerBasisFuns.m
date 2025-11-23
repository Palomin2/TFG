function ders = DerBasisFuns(i,x,p,n,knotVec) 
% ders(k,j) is the kth-1 derivative of the function 
% N_{i-p+j-1,p} where 1 <= k <= n+1 and 1 <= j <= p+1.
ders = zeros(n+1,p+1);
N = zeros(p+1,p+1);
N(1,1) = 1;
for j = 1:p
    left(j) = x - knotVec(i+1-j);
    right(j) = knotVec(i+j) - x;
    aux = 0;
    for r = 0:j-1
        N(j+1,r+1) = right(r+1)+left(j-r);
        temp = N(r+1,j)/N(j+1,r+1);
        N(r+1,j+1) = aux + right(r+1)*temp;
        aux = left(j-r)*temp;
    end
    N(j+1,j+1) = aux;
end
for j = 1:p+1
    ders(1,j) = N(j,p+1);
end
%
% Calculate the derivatives
%
for r = 1:p+1
    s1 = 0;
    s2 = 1;
    a(1,1) = 1;
    k = 1;
    d = 0;
    rk = r-1-k;
    pk = p-k;
    if r-1 >=k
        a(s2+1,1) = a(s1+1,1)/N(pk+2,rk+1);
        d = a(s2+1,1)*N(rk+1,pk+1);
    end
    if rk >=-1
        j1 = 1;
    else
        j1 = -rk;
    end
    if r-2 <=pk
        j2 = 0;
    else
        j2 = p-r-1;
    end
    for j=j1:j2
        a(s2+1,j+1) = (a(s1+1,j+1)-a(s1+1,j))/N(pk+2,rk+j+1);
        d = d + a(s2+1,j+1)*N(rk+j+1,pk+1);
    end
    if r-1 <= pk
        a(s2+1,k+1) = -a(s1+1,k)/N(pk+2,r);
        d = d + a(s2+1,k+1)*N(r,pk+1);
    end
    ders(k+1,r) = d;
    j = s1;
    s1 = s2;
    s2 = j;
end
r = p;
for j = 1:p+1
    ders(k+1,j) = ders(k+1,j)*r;
end

