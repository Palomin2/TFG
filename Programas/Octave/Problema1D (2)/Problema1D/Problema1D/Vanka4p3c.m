function p = Vanka4p3c(p,bp,nu,N,w,h)
%%
q=8;
for it = 1:nu
    for i = 4:3:N-3
        [Ap1, sm1]=Api(N,p,i-3,h);
        [Ap2, sm2]=Api(N,p,i-2,h);
        [Ap3, sm3]=Api(N,p,i-1,h);
        [Ap4, sm4]=Api(N,p,i,h);
        [Ap5, sm5]=Api(N,p,i+1,h);
        [Ap6, sm6]=Api(N,p,i+2,h);
        [Ap7, sm7]=Api(N,p,i+3,h);
        Ap=[Ap1;Ap2;Ap3;Ap4;Ap5;Ap6;Ap7];

        Ab = [sm1(q+1:q+7);sm2(q:q+6);sm3(q-1:q+5);sm4(q-2:q+4);sm5(q-3:q+3);sm6(q-4:q+2);sm7(q-5:q+1)]; 
        p([i-3:i+3]) = p([i-3:i+3]) + w*inv(Ab)*(bp([i-3:i+3])-Ap);
    end
    for i = 5:3:N-3
        [Ap1, sm1]=Api(N,p,i-3,h);
        [Ap2, sm2]=Api(N,p,i-2,h);
        [Ap3, sm3]=Api(N,p,i-1,h);
        [Ap4, sm4]=Api(N,p,i,h);
        [Ap5, sm5]=Api(N,p,i+1,h);
        [Ap6, sm6]=Api(N,p,i+2,h);
        [Ap7, sm7]=Api(N,p,i+3,h);
        Ap=[Ap1;Ap2;Ap3;Ap4;Ap5;Ap6;Ap7];
        
        Ab = [sm1(q+1:q+7);sm2(q:q+6);sm3(q-1:q+5);sm4(q-2:q+4);sm5(q-3:q+3);sm6(q-4:q+2);sm7(q-5:q+1)]; 
        p([i-3:i+3]) = p([i-3:i+3]) + w*inv(Ab)*(bp([i-3:i+3])-Ap);
    end
    for i = 6:3:N-3
        [Ap1, sm1]=Api(N,p,i-3,h);
        [Ap2, sm2]=Api(N,p,i-2,h);
        [Ap3, sm3]=Api(N,p,i-1,h);
        [Ap4, sm4]=Api(N,p,i,h);
        [Ap5, sm5]=Api(N,p,i+1,h);
        [Ap6, sm6]=Api(N,p,i+2,h);
        [Ap7, sm7]=Api(N,p,i+3,h);
        Ap=[Ap1;Ap2;Ap3;Ap4;Ap5;Ap6;Ap7];
        
        Ab = [sm1(q+1:q+7);sm2(q:q+6);sm3(q-1:q+5);sm4(q-2:q+4);sm5(q-3:q+3);sm6(q-4:q+2);sm7(q-5:q+1)]; 
        p([i-3:i+3]) = p([i-3:i+3]) + w*inv(Ab)*(bp([i-3:i+3])-Ap);
    end
end