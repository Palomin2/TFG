function p = Vanka3p2c(p,bp,nu,N,w,h)
%%
q=4;
for it = 1:nu
    for i = 3:2:N-2
        [Ap1, sm1]=Api(N,p,i-2,h);
        [Ap2, sm2]=Api(N,p,i-1,h);
        [Ap3, sm3]=Api(N,p,i,h);
        [Ap4, sm4]=Api(N,p,i+1,h);
        [Ap5, sm5]=Api(N,p,i+2,h);
        Ap=[Ap1;Ap2;Ap3;Ap4;Ap5;];
        
        Ab = [sm1(q+1:q+5);sm2(q:q+4);sm3(q-1:q+3);sm4(q-2:q+2);sm5(q-3:q+1)]; 
        p([i-2:i+2]) = p([i-2:i+2]) + w*inv(Ab)*(bp([i-2:i+2])-Ap);
    end
    for i = 4:2:N-2
        [Ap1, sm1]=Api(N,p,i-2,h);
        [Ap2, sm2]=Api(N,p,i-1,h);
        [Ap3, sm3]=Api(N,p,i,h);
        [Ap4, sm4]=Api(N,p,i+1,h);
        [Ap5, sm5]=Api(N,p,i+2,h);
        Ap=[Ap1;Ap2;Ap3;Ap4;Ap5;];
        
        Ab = [sm1(q+1:q+5);sm2(q:q+4);sm3(q-1:q+3);sm4(q-2:q+2);sm5(q-3:q+1)]; 
        p([i-2:i+2]) = p([i-2:i+2]) + w*inv(Ab)*(bp([i-2:i+2])-Ap);
    end
end