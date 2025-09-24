function p = Red_Black(A,p,bp,nu,N,w)

for it = 1:nu
    for j = 1:2:N
        p(j) = p(j) + w*(bp(j)-A(j,:)*p)/A(j,j);
    end
    for j = 2:2:N
        p(j) = p(j) + w*(bp(j)-A(j,:)*p)/A(j,j);
    end     
end