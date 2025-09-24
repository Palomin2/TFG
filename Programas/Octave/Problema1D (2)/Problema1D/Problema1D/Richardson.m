function p = Richardson(A,p,bp,nu,N,w)

w = 2/N*0.7; % Repasar w=0.7*2/(N-p);
for it = 1:nu
    pold = p;
    for i = 1:N
        p(i) = pold(i) + w*(bp(i)-A(i,:)*pold); % en lugar de dividir entre A(i,i), multiplicamos por un parametro w
    end
end