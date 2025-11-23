global q;
p = 2;
q = 2;

noPtsX = p+1;
noPtsY = q+1;

gcoord=meshRectangularCoord(2,2,noPtsX-1,noPtsY-1); 
controlPts=gcoord; % malla uniforme del cuadrado [0,2]
weights = ones(1,noPtsX*noPtsY)'; % pesos, tantos como length(controlPts)


knotUTemp = linspace(0,2,noPtsX-p+1); % Para completar el vector de nudos 
knotVTemp = linspace(0,2,noPtsY-q+1); % Para completar el vector de nudos

%%%% p=q=1
% uKnot = [0  knotUTemp 2]; 
% vKnot = [0  knotVTemp 2];

%% p=q=2
uKnot = [0 0 knotUTemp 2 2];
vKnot = [0 0 knotVTemp 2 2];

%%%% p=q=3
% uKnot = [0 0 0 knotUTemp 2 2 2];
% vKnot = [0 0 0 knotVTemp 2 2 2];

%%%% p=q=4
% uKnot = [0 0 0 0 knotUTemp 2 2 2 2]; 
% vKnot = [0 0 0 0 knotVTemp 2 2 2 2];

%%%% p=q=5
% uKnot = [0 0 0 0 0 knotUTemp 2 2 2 2 2]; 
% vKnot = [0 0 0 0 0 knotVTemp 2 2 2 2 2];

%%%% p=q=6
% uKnot = [0 0 0 0 0 0 knotUTemp 2 2 2 2 2 2]; 
% vKnot = [0 0 0 0 0 0 knotVTemp 2 2 2 2 2 2];

%%%% p=q=7
% uKnot = [0 0 0 0 0 0 0 knotUTemp 2 2 2 2 2 2 2]; 
% vKnot = [0 0 0 0 0 0 0 knotVTemp 2 2 2 2 2 2 2];

% %%% p=q=8
% uKnot = [0 0 0 0 0 0 0 0 knotUTemp 2 2 2 2 2 2 2 2]; 
% vKnot = [0 0 0 0 0 0 0 0 knotVTemp 2 2 2 2 2 2 2 2];
% noGPs=p+1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASO ANNULUS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % CASO ANNULUS
% % % clear controlPts
% % % clear weights
% % % 
% % % a = 0.3; % inner radius
% % % b = 0.5; % outer radius
% % % 
% % % % quadratic NURBS
% % % 
% % % p = 2;
% % % q = 2;
% % % noGPs = p+1;
% % % 
% % % uKnot = [0 0 0 1 1 1];
% % % vKnot = [0 0 0 1 1 1];
% % % clear controlPts
% % % % controlPts=[a 0;a a;0 a;
% % % %             0.5*(a+b) 0; 0.5*(a+b) 0.5*(a+b);0 0.5*(a+b); 
% % % %             b 0; b b; 0 b]';
% % % controlPts=[a 0; 0.5*(a+b) 0        ; b 0; 
% % %             a a; 0.5*(a+b) 0.5*(a+b); b b;
% % %             0 a; 0 0.5*(a+b)        ; 0 b];
% % % % figure(1)
% % % % plot(controlPts(:,1),controlPts(:,2),'oblack')
% % % % axis equal
% % % % axis([0 0.5 0 0.5])
% % % % pause
% % % 
% % % noPtsX = length(uKnot)-p-1;
% % % noPtsY = length(vKnot)-q-1;
% % % 
% % % weights = ones(1,noPtsX*noPtsY)';
% % % weights(4)=1/sqrt(2);
% % % weights(5)=1/sqrt(2);
% % % weights(6)=1/sqrt(2);
% % % % 
% % % pref=0; % FUNCIONA PARA CUALQUIER GRADO!
% % % for i=1:pref  % Sin embargo... es costoso hacerlo de esta manera, sobre todo para orden alto.
% % %     [p,noPtsX,uKnot,q,noPtsY,vKnot,controlPts,weights]=prefinement(p,noPtsX,uKnot,q,noPtsY,vKnot,controlPts,weights);
% % % end
% % % noCtrPts       = noPtsX * noPtsY;
% % % noDofs         = noCtrPts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


refineCount = 4; % como maximo haremos 8 ref en 2D
max_level = 2;

hRefinement2d
noGPs = p+1;  

% figure(2)
% plot(controlPts(:,1),controlPts(:,2),'oblack')
% axis equal
% axis([0 0.5 0 0.5])
% pause

generateIGA2DMesh;

noCtrPts = size(controlPts,1);% no of control points
%noElems  = size(elConn,1);    % no of elements
Nu = zeros(max_level+1,1);
Nu(1) = size(controlPts,1);

% initialization
K = sparse(noCtrPts,noCtrPts); % global stiffness matrix 
u = zeros(noCtrPts,1);        % displacement vector
f = zeros(noCtrPts,1);        % external force vector    

% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); %  quadrature
% Assembling system of equation
% Stiffness matrix and external force vector
% Loop over elements (knot spans)
for e=1:noElems  % VOY A IR ELEMENTO A ELEMENTO!!!
   idu    = index(e,1); % Accedo al indice del uspan de mi elemento
   idv    = index(e,2); % Accedo al indice del vspan de mi elemento
   xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
   etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]  
   
   sctr   = element(e,:); %  element scatter vector % indices de las funciones base no nulas
   nn     = length(sctr);
   pts    = controlPts(sctr,:); %ptos de control donde dichas funciones valen 1 resp.
   
   % loop over Gauss points 

    for gp=1:size(W,1)                        
      pt      = Q(gp,:);  % filas tienen las dos coordenadas del punto                         
      wt      = W(gp);                            
      
      % compute coords in parameter space
      Xi      = parent2ParametricSpace(xiE,pt(1)); %a partir del pto gaussiano hallo el parametrico
      Eta     = parent2ParametricSpace(etaE,pt(2)); %a partir del pto gaussiano hallo el parametrico
      J2      = jacobianPaPaMapping(xiE,etaE); % dxi / d\tilde{xi}
      
      %programa en c:
      %[dRdxi, dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights');
      
      %adaptacion a matlab: 
      [dRdxi, dRdeta] = NURBS2Ddersmatlab(Xi,Eta,p,q,uKnot,vKnot,weights');
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      % [dRdxi' dRdeta'] es una matriz (p+1)*(q+1) x 2
      jacob = pts'*[dRdxi' dRdeta'];  % dx / dxi
      J1    = det(jacob); % | dx / dxi |
      
      %%%%%%%%%%%%%%%%%%% Mi ejemplo era u(x,y)=xy(x-2)(y-2);
      i=findspan(noPtsX-1,p,Xi,uKnot);  
      j=findspan(noPtsY-1,q,Eta,vKnot); 
      N = BasisFuns(i+1,Xi,p,uKnot);
      M = BasisFuns(j+1,Eta,q,vKnot);
      AUX= kron(M,N).*weights(sctr(:))'/sum(kron(M,N).*weights(sctr(:))'); %/w
      X= AUX*controlPts(sctr,1:2);
      % u(x,y)=xy(x-2)(y-2)
      % bx=-2*X(1)*(X(1)-2)-2*X(2)*(X(2)-2); % u(x,y)=xy(x-2)(y-2)
      %bx=2*pi^2*sin(pi*X(1))*sin(pi*X(2));
      %bx=0;
      
      
      c=1;
      bx= 2*c^2*pi*pi*sin(c*pi*X(1))*sin(c*pi*X(2));
      
%       %%% ANNULUS: 
%       a=0.3;
%       b=0.5;
%       % PARA ANNULUS CON CONDICIONES DE FRONTERA DIRICHLET HOMOGENEAS
%       bx= 2*c^2*pi*pi*sin(c*pi*X(1))*sin(c*pi*X(2))*(X(1)^2+X(2)^2-a^2)*(X(1)^2+X(2)^2-b^2)-...
%           2*c*pi*cos(c*pi*X(1))*sin(c*pi*X(2))*(4*X(1)^3+4*X(1)*X(2)^2-2*X(1)*(a^2+b^2))+-...
%           2*c*pi*sin(c*pi*X(1))*cos(c*pi*X(2))*(4*X(2)^3+4*X(2)*X(1)^2-2*X(2)*(a^2+b^2))-...
%           4*sin(c*pi*X(1))*sin(c*pi*X(2))*(4*(X(1)^2+X(2)^2)-(a^2+b^2));
%       bx=0;
      f(sctr)=f(sctr)+ bx * AUX' * J1 * J2 *wt;
      %clear AUX
      %%%%%%%%%%%%%%%%%%%%%
      
      % Jacobian inverse and spatial derivatives
            
      dRdx       = [dRdxi' dRdeta']/jacob; %multiplico por la inversa del jacobiano 
                                           % es decir, divido entre dx/dxi
      % Asi obtenemos dR/dx = dR/dxi * (dx/dxi)^{-1} = dR/dx
      B       = dRdx;
      K(sctr,sctr) = K(sctr,sctr) + B * B' * J1 * J2 * wt;% gradiente por gradiente
                                                          % luego B*B'
      e
    end
end


% Solve the equation
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
                    % used to keep the  conditioning of the K matrix
%udofs  = [1 noCtrPts];  % global indecies  of the fixed x displacements

udofs  = [1:noPtsX];
%uFixed= [sin(pi*controlPts(1:noPtsX,1)).*sin(pi*controlPts(1:noPtsX,2))]';
uFixed=zeros(1,length(1:noPtsX));
for i=2:noPtsY-1
    udofs  = [udofs (i-1)*noPtsX+1 i*noPtsX];  % global indecies  of the fixed x displacements
    uFixed= [uFixed zeros(1,2)];
end 
udofs  = [udofs (((noPtsY-1)*noPtsX+1):noPtsY*noPtsX)];
uFixed  = [uFixed zeros(1,length((((noPtsY-1)*noPtsX+1):noPtsY*noPtsX)))]';
%uFixed=zeros(length(udofs),1); %%% CONDICION DE CONTORNO DIRICHLET HOMOGENEA

f=f-K(:,udofs)*uFixed;  % modify the  force vector

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(:,udofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
f(udofs)=bcwt*uFixed;

% SOLVE SYSTEM
U=K\f; % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PINTAR GRÁFICAS

% red=zeros(sqrt(length(U)),sqrt(length(U)));
% pesos=zeros(sqrt(length(U)),sqrt(length(U)));
% for i=1:sqrt(length(U))
%     red(i,:)=U(sqrt(length(U))*(i-1)+1:sqrt(length(U))*i);
%     pesos(i,:)=weights(sqrt(length(U))*(i-1)+1:sqrt(length(U))*i);
% end
% pesos
% pause
% x=zeros(size(red));
% y=zeros(size(red));
% 
% %for j=1:size(red,2)
%     for i=1:size(red,1)
% %         x(i,j)=2*(i-1)/(size(red,1)-1);
%         x(i,:)=controlPts((i-1)*(size(red,1))+1:i*(size(red,1)),1);
%     end
% %end
% %for i=1:size(red,1)
%     for j=1:size(red,2)
% %         y(i,j)=2*(j-1)/(size(red,2)-1);
% %         y(:,j)=controlPts(1:size(red):end,2);
%         y(j,:)=controlPts((j-1)*(size(red,1))+1:j*(size(red,1)),2);
%     end
% %end
% 
% P=zeros(size(red,1),size(red,2),4);
% P(:,:,1)=x;
% P(:,:,2)=y;
% P(:,:,3)=red;
% P(:,:,4)=pesos;
% Pw=P;
% xx=zeros(51,51);
% yy=zeros(51,51);
% zz=zeros(51,51);
% ww=zeros(51,51);
% for i=1:51
%     %%xx(i)=2*(i-1)/50;
%   for j=1:51
%       %%yy(j)=2*(j-1)/50;
%       %S=NURBSSurfacepoint(noPtsX,p,uKnot,noPtsY,q,vKnot,P,2*(i-1)/50,2*(j-1)/50);
%       S=NURBSSurfacepoint(noPtsX,p,uKnot,noPtsY,q,vKnot,P,2*(i-1)/50,2*(j-1)/50);
%       xx(i,j)=S(1);
%       yy(i,j)=S(2); 
%       %S=NURBSSurfacepoint(noPtsX,p,uKnot,noPtsY,q,vKnot,P,S(1),S(2));
%       zz(i,j)=S(3);
%       % ww(i,j)= xx(i)^2*yy(j)^2 -2*(xx(i)^2*yy(j)) -2*(xx(i)*yy(j)^2) +4*xx(i)*yy(j);
%       % ww(i,j)=sin(pi*(xx(i,j)))*sin(pi*(yy(i,j)));
%       %ww(i,j)=sin(pi*((i-1)/50))*sin(pi*((j-1)/50));
%       %ww(i,j)=xx(i,j);%*(2-xx(i,j));
%       ww(i,j)=sin(c*pi*xx(i,j))*sin(c*pi*yy(i,j))*(xx(i,j)^2+yy(i,j)^2-a^2)*(xx(i,j)^2+yy(i,j)^2-b^2);
%   end 
% end 
% figure(3)
% surface(xx,yy,zz);
%axis equal
%axis([0 0.5 0 0.5])

% figure(4)
% surface(xx,yy,zz-ww);

% figure(5)
% surface(xx,yy,ww);
% % axis equal
% % axis([0 0.5 0 0.5])

% error_inf=max(max(abs(zz-ww)))
% error_norma2=norm(zz-ww,2)

% figure(9)
% hold on
% plot(xx(1,1:51),yy(1,1:51),'black',xx(51,1:51),yy(51,1:51),'black',xx(1:51,1),yy(1:51,1),'black',xx(1:51,51),yy(1:51,51),'black')
% % PARA COMPARAR CUANDO HAGAMOS p-refinement
% % plot(savexx(1,1:51),saveyy(1,1:51),'blue',savexx(51,1:51),saveyy(51,1:51),'blue',savexx(1:51,1),saveyy(1:51,1),'blue',savexx(1:51,51),saveyy(1:51,51),'blue')
% hold off

