function p = Vanka4p4c(A,p,bp,nu,N,w)
%%%%
for it = 1:nu
    for i = 4:4:N-3
        Ab = A(i-3:i+3,i-3:i+3);
        p([i-3:i+3]) = p([i-3:i+3]) + w*inv(Ab)*(bp([i-3:i+3])-A([i-3:i+3],:)*p);
    end
    for i = 5:4:N-3
        Ab = A(i-3:i+3,i-3:i+3);
        p([i-3:i+3]) = p([i-3:i+3]) + w*inv(Ab)*(bp([i-3:i+3])-A([i-3:i+3],:)*p);
    end
    for i = 6:4:N-3
        Ab = A(i-3:i+3,i-3:i+3);
        p([i-3:i+3]) = p([i-3:i+3]) + w*inv(Ab)*(bp([i-3:i+3])-A([i-3:i+3],:)*p);
    end
    for i = 7:4:N-3
        Ab = A(i-3:i+3,i-3:i+3);
        p([i-3:i+3]) = p([i-3:i+3]) + w*inv(Ab)*(bp([i-3:i+3])-A([i-3:i+3],:)*p);
    end
end