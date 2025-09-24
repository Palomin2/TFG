function p = smooth(p,bp,nu,w,N,smth,h)
if smth == 1
    p = Gauss_Seidel(A,p,bp,nu,N,w);  
elseif smth == 2
    p = Red_Black(A,p,bp,nu,N,w); 
elseif smth == 3
    p = Jacobi(A,p,bp,nu,N,w); 
elseif smth == 4
    p = Richardson(A,p,bp,nu,N,w); 
elseif smth == 5
    p = Vanka2p(p,bp,nu,N,w,h); %funciona para todo p>=2
elseif smth == 6
    p = Vanka2p2colores(p,bp,nu,N,w,h); %funciona para todo p>=2
elseif smth == 7
    p = Vanka2p3colores(p,bp,nu,N,w,h); %funciona para todo p>=2
elseif smth == 8
    p = Vanka3p(p,bp,nu,N,w,h);    %funciona para todo p>=4
elseif smth == 9
    p = Vanka3p2c(p,bp,nu,N,w,h);  %funciona para todo p>=4
elseif smth == 10
    p = Vanka3p3c(p,bp,nu,N,w,h);  %funciona para todo p>=4
elseif smth == 11
    p = Vanka3p4c(p,bp,nu,N,w,h);   
elseif smth == 12
    p = Vanka4p(p,bp,nu,N,w,h);   %funciona para todo p>=6
elseif smth == 13
    p = Vanka4p2c(p,bp,nu,N,w,h);  %funciona para todo p>=6
elseif smth == 14
    p = Vanka4p3c(p,bp,nu,N,w,h);  %funciona para todo p>=6
elseif smth == 15
    p = Vanka4p4c(p,bp,nu,N,w,h);  
end
return