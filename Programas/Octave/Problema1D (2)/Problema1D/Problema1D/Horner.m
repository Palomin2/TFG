function C = Horner(a,n,u0)

%Evalua la curva C = \sum_{i=0}^{n} a_i * u^{i} en el punto u_0
% Los a[i] son las derivadas i-esimas de la curva en t=0 divididos entre i!
% Ejemplo: La curva (x-1)^2 tiene a=[1,-2,1]. Si la evaluamos en u0=1/2
% el resultado sera 1/4
C=a(n+1);
for i=n:-1:1
    C=C*u0 + a(i);
end 


%%%% PARA PINTAR UNA CURVA
% for i=1:101
%     x(i)=(i-1)/100;
% end 
% for i=1:101
%     y(i)=Horner(a,n,(i-1)/100);
% end 
% plot(x,y);