

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

p=2;
refinement=6;

controlPts = [  -1,  -1,  0, 1, 1;
                0,  1,  1,  1,  0]';

weights = [1, sqrt(2)/2, 1, sqrt(2)/2, 1];

knotVec    = [0, 0, 0, 1/2, 1/2, 1, 1, 1];

generateIGA1DMesh
clear weights
weights = controlPts(:,3);
noCtrPts = size(controlPts,1); % no of control points
noElems  = size(elConn,1);    % no of elements

% initialization
K = sparse(noCtrPts,noCtrPts); % global stiffness matrix
M = sparse(noCtrPts,noCtrPts); % matriz de masa
u = zeros(noCtrPts,1);        % displacement vector
f = zeros(noCtrPts,1);        % external force vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSamples  = 50000;
xi_samp   = linspace(0,1,nSamples);
Jvals     = zeros(size(xi_samp));
for k = 1:nSamples
    xi  = xi_samp(k);
    i   = findspan(noCtrPts-1,p,xi,knotVec);
    B   = BasisFuns(i+1,xi,p,knotVec);
    dB  = DerBasisFuns(i+1,xi,p,1,knotVec);
    conn = i-p+1:i+1;
    W    = sum(B.*weights(conn)');
    Wd   = sum(dB(2,:).*weights(conn)');
    dRdxi= ((dB(2,:).*weights(conn)')*W - Wd*(B.*weights(conn)'))/W^2;
    rxi  = (controlPts(conn,1:2))' * dRdxi';
    Jvals(k) = norm(rxi);
end
s_samp = cumtrapz(xi_samp,Jvals);     % longitud acumulada
L = s_samp(end);
s_samp = s_samp / L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noGPs=p+1;
% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 1 ); %  quadrature
% Assembling system of equation
% Stiffness matrix and external force vector
% Loop over elements (knot spans)
for e=1:noElems
   xiE   = elRange(e,:); % [xi_i,xi_i+1]
   conn  = elConn(e,:);
   noFns = length(conn);
   aaaaaaaaaaaa=111111
   % loop over Gauss points
    for gp=1:size(W,1)
      pt      = Q(gp,:);
      wt      = W(gp);
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * pt + xiE(2) + xiE(1)) % coord in parameter space
      J2      = 0.5 * ( xiE(2) - xiE(1) );

      i = findspan(noCtrPts-1,p,Xi,knotVec);
      wi=weights(conn);
      B=BasisFuns(i+1,Xi,p,knotVec);
      D = DerBasisFuns(i+1,Xi,p,2,knotVec)
      B1=D(2,:);
      B2=D(3,:);

      w=0;
      w=sum(BasisFuns(i+1,Xi,p,knotVec).*wi');
      N = B.*wi'/w;
      W1= sum(B1.*wi');
      dNdxi=((B1.*wi')*w - W1*B.*wi')/w^2;

      dN2    = DerBasisFuns(i+1,Xi,p,2,knotVec); % 2ª derivada B-spline
      W2 = sum(B2.*wi');                         % W''

      dN2dxi2 = ((B2.*wi')*w^2 - 2*W1*w*(B1.*wi') - W2*w*(B.*wi') + 2*W1^2*(B.*wi')) / w^3;

      jacob1 = dNdxi*controlPts(conn,1:2)
      J1     = norm(jacob1)
      dNdx   = (1/J1)*dNdxi;

      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      K(conn,conn) = K(conn,conn) + dNdx' * dNdx * J1 * J2 * wt;

      % compute the external force, kind of body force ???
      s_gp = interp1(xi_samp, s_samp, Xi, 'cubic');
      c=10;
      f_gp = (c*pi)^2 * sin(c*pi * s_gp)/pi^2; % DIVIDO ENTRE PI^2 POR EL REESCALADO DE LA LONGITUD PI DE LA CURVA
      wt
      J1
      J2
      N'
      N' * f_gp * J1 * J2 * wt
      f(conn) = f(conn) + N' * f_gp * J1 * J2 * wt;
    end
end
%
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

% SOLVE SYSTEM

U=K\f;

a=0:1/500:1;
ordenadas2 = sin(c*pi*a);

xxx=-cos(pi*a);
yyy=sin(pi*a);
zzz=sin(c*pi*a);

figure(2)
aux2=zeros(500,4);
for j=1:501
    s_gp = interp1(xi_samp, s_samp,(j-1)/500, 'cubic');

    C=NURBSCurvepoint3D(length(controlPts),p,knotVec,[controlPts(:,1),controlPts(:,2),U,controlPts(:,3)],(j-1)/500);
    % % C=dersNURBSCurvepoint3D(length(controlPts),p,knotVec,[controlPts(:,1),controlPts(:,2),U,controlPts(:,3)],(j-1)/500);

    zzz(j)=sin(c*pi*s_gp);
    % zzz(j)=c*pi*cos(c*pi*s_gp);

    xx(j)=C(1);
    yy(j)=C(2);
    zz(j)=C(3);
end

plot3(controlPts(:,1),controlPts(:,2),U,'black',xx,yy,zz,'blue',xx,yy,zzz,'red')


% plot3(controlPts(:,1),controlPts(:,2),U,'black',xx,yy,ordenadas,'red',xx,yy,zz,'blue',a,sqrt(1-a.^2),ordenadas2,'yellow')

Error_curve1D

return
