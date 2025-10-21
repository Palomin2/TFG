function [Ubar,Qw] = RefineKnotVectCurve(n,p,U,Pw,X,r)
%--------------------------------------------------------------
%function [Ubar,Qw] = RefineKnotVectCurve(n,p,U,Pw,X,r)
% NURBS-Book (algorithm A5.4) (modified)
% insert multiple knots into curve
%INPUT:
% n         : number ob basis functions -1 !
%        NURBS-Book: n+1 # basis, np max index (startindex 0)
%        here        n   # basis and max index (startindex 1)
% p          : degree of the basis functions
% U         : old knotvector
% Pw         : old control points
% X          : vector of new knots (multiple entries possible)
% r          :  (size of X) -1 (count the multple entries as well
%             reason: same as above: max X index
%OUTPUT:
% Ubar        : newknot vector
% Qw         : new control points
%--------------------------------------------------------------

%initialise arrays;
dim  = size(Pw,2);
Qw   = zeros(n+r+2,dim);
Ubar = zeros(1,n+p+1+r);
% n ES EL NUMERO DE FUNCIONES BASE -1
m = n+p+1;
a = findspan(n,p,X(1),U);
b = findspan(n,p,X(r+1),U);
b = b+1; % pues b es el indice del nudo u_b que es > que todos los x_i
% En los comentarios los indices empiezan en 0
for j=0:a-p % en 1 knot insertion ya ves que los alfa_i con i<=a-p son 1
    Qw(j+1,:) = Pw(j+1,:); % Q_{0},...,Q_{a-p} = P_{0},...,P_{a-p}
end
% Notar que los ptos desde Q_{b,r} a Q_{b+r-1,r} se modifican debido a que en
% la recurrencia de (5.15) el segundo sumando es Q_{i-1,r-1}.
% Por ejemplo: Q_{b+1,r} cambiará porque depende de Q_{b,r-1}.
for j=b-1:n % de Q_{a-p} a Q_{b+r-1} son los ptos de control modificados
    Qw(j+r+2,:) = Pw(j+1,:); % Q_{b+r},...,Q_{n+r+1} = P_{b-1},...,P_{n}
end
for j=0:a
    Ubar(j+1)= U(j+1); % ubar(0:a)= u(0:a)
end
for j=b+p:m % u(b+p) es el primer nudo que ya no cambia, b es el knot span  
    Ubar(j+r+2) = U(j+1); % ubar(b+p+r+1) = u(b+p), nuestro X es de longitud r+1
end
i=b+p-1;
k=b+p+r;
for j=r:-1:0
    while (X(j+1) <= U(i+1) && i>a) % Caso en el que quedan nudos antiguos que superen el actual nudo nuevo
        Qw(k-p,:) = Pw(i-p,:); % Por tanto se asigna por el final un pto de control antiguo 
        Ubar(k+1) = U(i+1); % Lo mismo con un nudo
        k=k-1; % Pues hay un nuevo nudo y pto de control asignados
        i=i-1; % Compararemos con el siguiente (anterior) nudo antiguo
    end
    % Ahora X(j+1) > U(i+1), luego X(j+1) será el siguiente nudo a
    % insertar, lo que supone un nuevo pto de control construido con los
    % alfas.
    
    % para k-p, alfa=1 y el pto de control va a ser "el posterior" porque
    % no cambia:
    Qw(k-p,:) = Qw(k-p+1,:); % El nuevo nudo va a a ser el menor valor encontrado en el espacio
                             % paramétrico, luego nuestro antiguo primer
                             % pto de control de supuesta multiplicidad p+1
                             % ha de seguir intacto. Para ello traslado los
                             % indices provocando un cambio de indices para
                             % expresar la formula 5.15.
    % Desde k-p+1 hasta k los alfas no son 0 ni 1, y los ptos de control
    % cambian. Notar que todas estas variables las estamos modificando
    % debido al efecto que tiene la insertación de un nudo nuevo, muchos de
    % los ptos de control que vamos a modificar ahora se modifican de nuevo
    % si en la siguiente iteración metemos otro nudo nuevo.
    for l=1:p % En este bucle calculamos los p puntos de control nuevos tras un knot insertion
        ind = k-p+l;
        alfa = Ubar(k+l+1) - X(j+1); % numerador de 1-alfa, pues vamos a añadir X(j+1) como nuevo nudo
        if (abs(alfa) == 0)
            Qw(ind,:) = Qw(ind+1,:); % Trivialidad... formula 5.15 como siempre.
        else
            alfa = alfa/(Ubar(k+l+1) - U(i-p+l+1)); % 1-alfa construido
            Qw(ind,:) = alfa* Qw(ind,:) + (1-alfa)* Qw(ind+1,:); % Es 5.15 pero cambiado el orden de suma
        end
    end
    Ubar(k+1) = X(j+1); % pues hemos superado el bucle while
    k=k-1; % Pues hay un nuevo nudo y pto de control asignados
end
end