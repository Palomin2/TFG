function N = BasisFuns(i,x,p,knotVec) 

% To calculate N_{i-p,p}(x), ...N_{i,p}(x)
% N(j) corresponds to N_{i-p+j-1}, j = 1,...,p+1.


N(1) = 1;
for j = 1:p
    left(j) = x - knotVec(i+1-j);
    right(j) = knotVec(i+j) - x;
    aux = 0;
    for r = 0:j-1
        temp = N(r+1)/(right(r+1)+left(j-r));
        N(r+1) = aux + right(r+1)*temp;
        aux = left(j-r)*temp;
    end
    N(j+1) = aux;
end
