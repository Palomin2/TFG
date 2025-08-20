function N = BasisFuns(i,x,p,knotVec) 

% To calculate N_{i-p,p}(x), ...N_{i,p}(x)
% N(j) corresponds to N_{i-p+j-1}, j = 1,...,p+1.


N(1) = 1; % N_{i,0}
% En cada iteracion de j se hallan los j+1 B-splines no nulos de grado j en
% [t_j, t_{j+1}).
% Notar que en cada j-iteración se modifican por tanto los N(:).
% Al acabar cada iteración obtenemos:
% N_{i-j,j}, N_{i-j+1,j}, ... , N_{i,j}
% Hasta que llegamos a la iteración j=p.
for j = 1:p
    left(j) = x - knotVec(i+1-j);
    right(j) = knotVec(i+j) - x;
    aux = 0;
    for r = 0:j-1 % Calculamos N_{i-j+r,j}. Se empieza obteniendo N_{i-j,j} y acabamos hallando N_{i-1,j}
        temp = N(r+1)/(right(r+1)+left(j-r)); 
        N(r+1) = aux + right(r+1)*temp;
        aux = left(j-r)*temp; % Esta es la funcion base (multiplicada por el cociente) que no cuenta para calcular
                              % la N(r+1) en esta iteración, pero es
                              % necesaria en la siguiente por la relacion de recurrencia, por eso la
                              % guardamos.
    end
    N(j+1) = aux; % N_{i,j}
end
