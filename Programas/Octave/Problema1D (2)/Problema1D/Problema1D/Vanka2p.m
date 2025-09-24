function p = Vanka2p(p,bp,nu,N,w,h)
% global q; % PUEDE QUE PONER AQUI global q no sea eficiente, poner q=2 en ese caso
q=2;
for it = 1:nu
    for i = 2:N-1
        
        [Ap1, sm1]=Api(N,p,i-1,h);
        [Ap2, sm2]=Api(N,p,i,h);
        [Ap3, sm3]=Api(N,p,i+1,h);  %%%  %Ab = [sm1(4:6);sm2(3:5);sm3(2:4)];  q=3
        Ap=[Ap1;Ap2;Ap3];
        
        Ab = [sm1(q+1:q+3);sm2(q:q+2);sm3(q-1:q+1)]; %Ab = [sm1(3:5);sm2(2:4);sm3(1:3)];  q=2
        p([i-1:i+1]) = p([i-1:i+1]) + w*inv(Ab)*(bp([i-1:i+1])-Ap);
    end
end

