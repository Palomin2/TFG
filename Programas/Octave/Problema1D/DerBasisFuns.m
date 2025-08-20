function ders = DerBasisFuns(i,x,p,n,knotVec) 
% ders(k,j) is the kth-1 derivative of the function 
% N_{i-p+j-1,p} where 1 <= k <= n+1 and 1 <= j <= p+1.

% Realmente este programa solo calcula la 1a derivada
ders = zeros(n+1,p+1);
N = zeros(p+1,p+1);
N(1,1) = 1;
%Calculamos las funciones base no nulas en nuestro knot span de x
for j = 1:p
    left(j) = x - knotVec(i+1-j);
    right(j) = knotVec(i+j) - x;
    aux = 0;
    for r = 0:j-1
        N(j+1,r+1) = right(r+1)+left(j-r); % triangular inferior
        temp = N(r+1,j)/N(j+1,r+1);
        N(r+1,j+1) = aux + right(r+1)*temp; % triangular superior
        aux = left(j-r)*temp;
    end
    N(j+1,j+1) = aux;
end
for j = 1:p+1
    ders(1,j) = N(j,p+1); % la primera fila de ders son las funciones sin derivar
end
%
% Calculate the derivatives
%

%k=1, calculamos la primera derivada de todas las funciones en ders(2,:) 
%base
for r = 1:p+1
    s1 = 0;
    s2 = 1;
    a(1,1) = 1;
    k = 1;
    d = 0;
    rk = r-1-k;
    pk = p-k;
    
    %%% PRIMER SUMANDO EN LA EXPRESION DE N_{i-p+r,p}' :
    
    if r-1 >=k % si rk >= 0. Esto se cumple a partir de N_{i-p+1, p}' hasta N_{i, p}'
        a(s2+1,1) = a(s1+1,1)/N(pk+2,rk+1); % En N(pk+2,rk+1) se guarda u_{i+r-k+1}-u_{i-p+r}
        d = a(s2+1,1)*N(rk+1,pk+1); % En N(rk+1,pk+1) se guarda N_{i-p+r,p-k}
    end
    %%%
    if rk >=-1 % Para k=1, se cumple siempre.
        j1 = 1; 
    else
        j1 = -rk; % j1>1, pero no va a pasar
    end
    if r-2 <=pk % para k=1, esto ocurre siempre....
        j2 = 0;
    else   % este else nos sobra, nunca va a pasar con k=1
        j2 = p-r-1; 
    end
    
    %%%SUMANDOS INTERMEDIOS:
    %%% Se necesitarían para derivadas de mayor orden.
    for j=j1:j2 %como j2=0 siempre con k=1, este for nos sobra
        a(s2+1,j+1) = (a(s1+1,j+1)-a(s1+1,j))/N(pk+2,rk+j+1);
        d = d + a(s2+1,j+1)*N(rk+j+1,pk+1);
    end
    %%%
    
    %%% ÚLTIMO SUMANDO EN LA EXPRESION DE N_{i-p+r,p}' :
    %%% el ultimo es N_{i-p+r+k-1, p-k} (debido a que r empieza tomando valor 1...)
    %%% La condicion para que esta función base sea no nula en [u_{i+1},u_i)
    %%% es que (i-p+r+k-1)<= i   ==> r-1 <= p-k.
    if r-1 <= pk % esto pasa siempre excepto con r=p+1, pues k=1
        a(s2+1,k+1) = -a(s1+1,k)/N(pk+2,r); % calculamos los "a_{k,k}"
        d = d + a(s2+1,k+1)*N(r,pk+1); % En N(r,pk+1) se guarda N_{i-p+r+k-1,p-k}
    end
    %%% 
    
    ders(k+1,r) = d; % directamente podiamos poner ders(2,r)=d
    
    j = s1; % estos intercambios de filas en a[][] nos sobran 
    s1 = s2; % los necesitariamos si hubiera un bucle en k
    s2 = j; % pero solo nos preocupamos por k=1, luego se pueden quitar.
end
r = p;
for j = 1:p+1
    ders(k+1,j) = ders(k+1,j)*r; % esto nos sirve solo con k=1, sino mirar NURBS book: p!/(p-k)!
end

