function p = Jacobi(A,p,bp,nu,N,w)

for it = 1:nu
    pold = p;
    for i = 1:N
        p(i) = pold(i) + w*(bp(i)-A(i,:)*pold)/A(i,i);
    end
end