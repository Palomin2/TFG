p = 4 ; %%%%%% cambiar en el programa Api y en tu suavizador!!!!!!!!!
global q;
q=p;
% Linear
% knotVec    = [0 0 2 2];
% controlPts = [0 0; 2 0];
% p = 1;
% noGPs = 2;
% --------------------------
if p == 2
% Quadratic
    knotVec    = [0 0 0 2 2 2];
    controlPts = [0 0;1 0; 2 0];
    noGPs   = 3;

% --------------------------
% Cubic
elseif p == 3
    knotVec    = [0 0 0 0 2 2 2 2];
    controlPts = [0 0;2/3 0;4/3 0;2 0];
    noGPs = 4;
% --------------------------
%Cuarto orden
elseif p == 4
    knotVec    = [0 0 0 0 0 2 2 2 2 2];
    controlPts = [0 0;0.5 0;1 0;1.5 0; 2 0];
    noGPs = 5;
% --------------------------
%Quinto orden
elseif p == 5
    knotVec    = [0 0 0 0 0 0 2 2 2 2 2 2];
    controlPts = [0 0;2/5 0;4/5 0;6/5 0;8/5 0; 2 0];
    noGPs = 6;
% --------------------------
%Sexto orden
elseif p == 6
    knotVec    = [0 0 0 0 0 0 0 2 2 2 2 2 2 2];
    controlPts = [0 0;2/6 0;2/3 0;1 0;4/3 0;10/6 0; 2 0];
    noGPs = 7;
% --------------------------
%Septimo orden
elseif p == 7
    knotVec    = [0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2];
    controlPts = [0 0;2/7 0;4/7 0;6/7 0;8/7 0;10/7 0;12/7 0;2 0];
    noGPs = 8;
% Octavo orden
elseif p == 8
    knotVec    = [0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2];
    %controlPts = [0 0;2/8 0;4/8 0;6/8 0;1 0;10/8 0;12/8 0;14/8 0;2 0];
    controlPts = [0 0;2/8 0;4/8 0;6/8 0;1 0;10/8 0;12/8 0;14/8 0;2 0];
    noGPs = 9;
end
%Noveno orden
% knotVec    = [0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2];
% controlPts = [0 0;2/9 0;4/9 0;6/9 0;8/9 0;10/9 0;12/9 0;14/9 0;16/9 0;2 0];
% p = 9;
% noGPs = 10;

weights = ones(1,size(controlPts,1));
% weights = [1,20,1];
% weights = [1,1,1,1,10,1,1,1,1];
refinement =5;
max_level = 9; % maximo numero de niveles en la jerarquia del multigrid

WI=cell(max_level,1);   %%%%% Para la restriccion con NURBS %%%%%
WI{1}=weights;          %%%%%                               %%%%%

generateIGA1DMesh
clear weights
weights = controlPts(:,3);
%controlPts
% pause
noCtrPts = size(controlPts,1); % no of control points
noElems  = size(elConn,1);    % no of elements
Nu = zeros(max_level+1,1); %numero de ptos de control en cada malla
Nu(1) = size(controlPts,1)
% pause
if p == 2
    for j = 2:max_level
        Nu(j) = Nu(j-1)/2+1;
    end
elseif p == 3
    for j = 2:max_level
        Nu(j) = (Nu(j-1)+1)/2+1;
    end
elseif p == 4
    for j = 2:max_level
        Nu(j) = Nu(j-1)/2+2;
    end
elseif p == 5
    for j = 2:max_level
        Nu(j) = (Nu(j-1)+1)/2+2;
    end
elseif p == 6
    for j = 2:max_level
        Nu(j) = Nu(j-1)/2+3;
    end
elseif p == 7
    for j = 2:max_level
        Nu(j) = (Nu(j-1)+1)/2+3;
    end
elseif p == 8
    for j = 2:max_level
        Nu(j) = Nu(j-1)/2+4;
    end
end
% initialization
K = sparse(noCtrPts,noCtrPts); % global stiffness matrix
M = sparse(noCtrPts,noCtrPts); % matriz de masa
u = zeros(noCtrPts,1);        % displacement vector
f = zeros(noCtrPts,1);        % external force vector
% pause
% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 1 ); %  quadrature
% Assembling system of equation
% Stiffness matrix and external force vector
% Loop over elements (knot spans)
for e=1:noElems
   xiE   = elRange(e,:); % [xi_i,xi_i+1]
   conn  = elConn(e,:);
   noFns = length(conn);

   % loop over Gauss points
    for gp=1:size(W,1)
      pt      = Q(gp,:);
      wt      = W(gp);
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * pt + xiE(2) + xiE(1)); % coord in parameter space
      J2      = 0.5 * ( xiE(2) - xiE(1) );

      i = findspan(noCtrPts-1,p,Xi,knotVec);
      %N = BasisFuns(i+1,Xi,p,knotVec);
      w=0;
      w=sum(BasisFuns(i+1,Xi,p,knotVec).*weights(conn)');
      N = BasisFuns(i+1,Xi,p,knotVec).*weights(conn)'/(w);
      dNdx = DerBasisFuns(i+1,Xi,p,1,knotVec);
      derw= sum(dNdx(2,:).*weights(conn)');
      dNdxi=((dNdx(2,:).*weights(conn)')*w - derw*(BasisFuns(i+1,Xi,p,knotVec).*weights(conn)'))/w^2;
      %dNdx = DerBasisFuns(i+1,Xi,p,1,knotVec);
      %dNdxi = dNdx(2,:);
      %jacobN=dNdx(1,:)*controlPts(conn,1:2); % MATRIZ DE MASA
      dMdx=N;                % MATRIZ DE MASA
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      jacob1 = dNdxi*controlPts(conn,1:2);
      J1     = norm (jacob1);
      dNdx   = (1/J1)*dNdxi;


      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      M(conn,conn) = M(conn,conn) + dMdx' * dMdx * J2 * wt;% MATRIZ DE MASA
      K(conn,conn) = K(conn,conn) + dNdx' * dNdx * J1 * J2 * wt;

      % compute the external force, kind of body force

      X       = N * controlPts(conn,1:2);
      % u(x) correspondiente al ejemplo de Hughes
%       bx      = X(1);
      c=10;
      bx      = c^2*pi^2*sin(c*pi*X(1));
      f(conn) = f(conn) + bx * N' * J1 * J2 * wt;
    end
end


% Solve the equation
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
                    % used to keep the  conditioning of the K matrix
udofs  = [1 noCtrPts];  % global indecies  of the fixed x displacements
uFixed = [0 0]'; % BCs: u[0]=u[1]=0
%uFixed = [0 1]'; % BCs: u[0]=0;u[1]=1;

f=f-K(:,udofs)*uFixed;  % modify the  force vector

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(:,udofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
f(udofs)=bcwt*uFixed;

%pause

% SOLVE SYSTEM

U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PINTAR SOLUCION

aux=zeros(500,3);
for j=1:501
    C=NURBSCurvepoint(length(U),p,knotVec,[controlPts(:,1),U,controlPts(:,3)],2*(j-1)/500);
    aux(j,:)=C;
    xx(j)=C(1);
    yy(j)=C(2);
end

plot(xx,yy);
figure(1)
%  abcisas = 0:1/200:1 ;
%  ordenadas = (-abcisas.^3 + abcisas)/6 ;
ordenadas = sin(c*pi*xx) ;
% ordenadas=xx.*(2-xx);
plot(controlPts(:,1),U,'black',xx,ordenadas,'red',xx,yy,'blue')

% max(abs(yy-ordenadas))
% norm(yy-ordenadas,2)

Error1D
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear w
%  Multigrid parameters
%

tol = 10^(-10);
smth=10; % tipo de suavizador % 7 10 14 son los de 3colores
w = 1; % parametro de relajacion (si hace falta)
nu1 = 1; % numero de presuavizados
nu2 = 0; % numero de postsuavizados
cyc_ind = 1;  % 1 for V-cycle, and 2 for W-cycle

R = cell(max_level,1); % matrices de restriccion
A = cell(max_level,1); % matrices de rigidez
P = cell(max_level,1); % matrices de prolongacion
B = cell(max_level,1);

A{1} = K;
% B{1} = K;
for level = 2:max_level
    Nf = Nu(level-1); % numero de puntos de control del nivel anterior
    Nc = Nu(level); % numero de puntos de control del nivel actual
%     R{level} = restriction(Nf,p,Nc); % Calcula la matriz de restriccion
    R{level} = diag(WI{length(WI)-level+1})*restriction(Nf,p,Nc)*inv(diag((WI{length(WI)-level+2})));
%     R1 = R{level};
%     full(R1)
%     pause
    P{level} = R{level}'; % la matriz de prolongacion es la traspuesta de la restriccion
%     A{level} = R{level}*A{level-1}*P{level}; % matriz de rigidez de la malla level a partir
                                             % de la de level-1, que es m s
                                             % fina. Se calcula con:
                                             % Restriccion*Rigidez*Prolongacion
%     A{level}=Matrices(p,refinement-level+1);
end

res = 1.0; % inicializacion del residuo
res_old = 1.0; % inicializacion del residuo
u = rand(noCtrPts,1); % aproximacion inicial del sistema
% u = 10*ones(noCtrPts,1); % aproximacion inicial del sistema
res_in_p = max(abs(f-K*u)); % residuo inicial
iter = 1;
h=1/2^(refinement-1);
% f = zeros(size(f));
his=[];
tic
while res>tol*res_in_p

% for iter = 1:70
      u = MG(R,P,u,f,nu1,nu2,Nu,w,cyc_ind,smth,max_level,1,h);  % Programa principal de multigrid
%       rp = f - K*u;
      for i=1:Nu(1)
          rp(i)=f(i)-Api(Nu(1),u,i,h); % HACER ESTO EMPEORA EL PROGRAMA?
      end
%       plot(rp)
      res = max(abs(rp)) % norma infinito del residuo
      iter=iter+1;
       res/res_old % cociente entre el residuo actual y el anterior para aproximar el factor
                   % de convergencia asintotico
      factor(iter) = res/res_old;
      res_old = res;
      his=[his factor(iter)];
      iter
%       plot(rp)
%       pause
end
toc
% pause
% max(abs(U-u)); % comparamos la solucion del sistema con la aproximacion de MG
% rho_m = mean(factor(40:70)) %hacemos la media de los ultimos factores para aproximar el factor de
                            % de convergencia


