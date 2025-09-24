function p = Vanka3p(A,p,bp,nu,N,w)

for it = 1:nu
    for i = 3:4:N-2
        Ab = A(i-2:i+2,i-2:i+2);
        p([i-2:i+2]) = p([i-2:i+2]) + w*inv(Ab)*(bp([i-2:i+2])-A([i-2:i+2],:)*p);
    end
    for i = 4:4:N-2
        Ab = A(i-2:i+2,i-2:i+2);
        p([i-2:i+2]) = p([i-2:i+2]) + w*inv(Ab)*(bp([i-2:i+2])-A([i-2:i+2],:)*p);
    end
    for i = 5:4:N-2
        Ab = A(i-2:i+2,i-2:i+2);
        p([i-2:i+2]) = p([i-2:i+2]) + w*inv(Ab)*(bp([i-2:i+2])-A([i-2:i+2],:)*p);
    end
    for i = 6:4:N-2
        Ab = A(i-2:i+2,i-2:i+2);
        p([i-2:i+2]) = p([i-2:i+2]) + w*inv(Ab)*(bp([i-2:i+2])-A([i-2:i+2],:)*p);
    end
end