function p = MG(R,P,p,bp,nu1,nu2,N,w,cyc_ind,smth,max_level,level,h)

% A: todas las matrices de rigidez correspondientes a los niveles de la
% jerarquia
% R: todas las matrices de restriccion correspondientes a los niveles de la
% jerarquia
% P: todas las matrices de prolongacion correspondientes a los niveles de la
% jerarquia
% p= Como argumento es la aproximacion en la iteracion actual, el programa
% devuelve la siguiente aproximacion.
% bp= rhs del nivel actual
% nu1= numero de presuavizados
% nu2= numero de postsuavizados
% N = numero de ptos de control en cada nivel, y por tanto tamaño del
% sistema en la iteracion actual
% w=parametro de relajacion
% cyc_ind= V ciclo o W ciclo
% smth = suavizador: GS, Jacobi...
% max_level= numero de niveles en la jerarquia
% level = nivel actual

% Performs a multigrid cycle; recursive.
% Af = A{level};
if (level == max_level)
%     p = Af\bp;
    p = smooth(p,bp,100,w,N(level),smth,h); % PONER ESTA PARTE
    return
end

Rf = R{level+1}; % Matriz de Restriccion al siguiente nivel
p = smooth(p,bp,nu1,w,N(level),smth,h); % suavizado de tu aproximacion
%%%%%%%%%%CALCULO DEL RESIDUO%%%%%%%%%%%%%

for i=1:N(level)
    rp(i)= bp(i) - Api(N(level),p,i,h);
end

% rp = bp - Af*p; % calculo del residuo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ncoarse = N(level+1); % numero de ptos de control en el siguiente nivel
bcp = Rf*rp'; % restriccion del residuo
% restricionrp=bcp(1);
% restricionrp
% pause  
pc = zeros(Ncoarse,1); % inicializamos a 0 la aproximacion en la siguiente malla
for rec_calls = 1:cyc_ind    
    pc = MG(R,P,pc,bcp,nu1,nu2,N,w,cyc_ind,smth,max_level,level+1,2*h);    % Recursive call
end
Pf = P{level+1}; %prolongacion al siguiente nivel (subimos de nivel)
p = p + Pf*pc; % corregir la aproximacion en el siguiente nivel sumandole la aproximacion en el nivel anterior
               % multiplicada por la prolongacion
p = smooth(p,bp,nu2,w,N(level),smth,h); % postsuavizado en cada nivel 

return
