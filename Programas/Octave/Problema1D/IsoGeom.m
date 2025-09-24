p = 4 ; %%%%%% cambiar en el programa Api!!!!!!!!!
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


generateIGA1DMesh
clear weights
weights = controlPts(:,3);
%controlPts
% pause
noCtrPts = size(controlPts,1); % no of control points
noElems  = size(elConn,1);    % no of elements

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
      c=1;
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

max(abs(yy-ordenadas))
norm(yy-ordenadas,2)
% % %
% % % pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

