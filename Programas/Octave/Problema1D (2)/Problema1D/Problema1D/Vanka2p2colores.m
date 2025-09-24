function p = Vanka2p2colores(p,bp,nu,N,w,h)
q=5;
for it = 1:nu
    for i = 2:2:N-1
        [Ap1, sm1]=Api(N,p,i-1,h);
        [Ap2, sm2]=Api(N,p,i,h);
        [Ap3, sm3]=Api(N,p,i+1,h);  
        Ap=[Ap1;Ap2;Ap3];
        
        Ab = [sm1(q+1:q+3);sm2(q:q+2);sm3(q-1:q+1)]; 
        p([i-1:i+1]) = p([i-1:i+1]) + w*inv(Ab)*(bp([i-1:i+1])-Ap);
    end
    for i = 3:2:N-1
        [Ap1, sm1]=Api(N,p,i-1,h);
        [Ap2, sm2]=Api(N,p,i,h);
        [Ap3, sm3]=Api(N,p,i+1,h);  
        Ap=[Ap1;Ap2;Ap3];
        
        Ab = [sm1(q+1:q+3);sm2(q:q+2);sm3(q-1:q+1)]; 
        p([i-1:i+1]) = p([i-1:i+1]) + w*inv(Ab)*(bp([i-1:i+1])-Ap);
    end
end

% En cada iteracion, modificamos las componentes i-1,i,i+1 de la
% aproximacion mediante un Gauss-Seidel por bloques 3x3. En Gauss-Seidel
% empleabamos el inverso del elemento diagonal, aquí tenemos que calcular
% la inversa de la matriz 3x3.

% De nuevo es la ecuacion del residuo por bloques, donde el residuo tiene 3
% componentes:  (bp([i-1:i+1])-A([i-1:i+1],:)*p)