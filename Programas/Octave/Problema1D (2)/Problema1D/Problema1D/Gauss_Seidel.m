function p = Gauss_Seidel(A,p,bp,nu,N,w)

for it = 1:nu
    for i = 1:N
        p(i) = p(i) + w*(bp(i)-A(i,:)*p)/A(i,i);
    end
end
    %Empleo la ecuacion del residuo
    %Notar que en p estamos empleando las aproximaciones anteriores, ya
    %estan guardadas.
            