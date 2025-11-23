% uuu=[];
% vvv=[];
a=2
b=0
L2_error=0;
L2_dererror=0;
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

%       xx=0;
%       yy=0;
%       zz=0;
%       Sderu=0;
%       Sderv=0;
%       exacta=0;
%       derexactau=0;
%       derexactav=0;
      S=NURBSSurfacepoint(noPtsX,p,uKnot,noPtsY,q,vKnot,Pw,Xi,Eta);
      xx=S(1);
      yy=S(2);
      zz=S(3);
      %Sder=dersNURBSSurfacepoint(noPtsX,p,uKnot,noPtsY,q,vKnot,Pw,Xi,Eta);
%       Sderu=Sder(2,3);
%       Sderv=Sder(1,3);
      %Sderu=Sder(1);
      %Sderv=Sder(2);
      exacta=sin(c*pi*xx)*sin(c*pi*yy)*(xx^2+yy^2-a^2)*(xx^2+yy^2-b^2);
      %derexactau=c*pi*cos(c*pi*xx)*sin(c*pi*yy)*(xx.^2*(xx.^2+2*yy.^2-b^2-a^2)+yy.^2*(yy .^2-b^2-a^2)+(a*b)^2) + ...
                       %2*sin(c*pi*xx)*sin(c*pi*yy)*(2*xx.^3+2*xx*yy.^2-xx*(a^2+b^2));
      %derexactav=c*pi*sin(c*pi*xx)*cos(c*pi*yy)*(xx.^2*(xx.^2+2*yy.^2-b^2-a^2)+yy.^2*(yy.^2-b^2-a^2)+(a*b)^2) + ...
                       %2*sin(c*pi*xx)*sin(c*pi*yy)*(2*yy.^3+2*yy*xx.^2-yy*(a^2+b^2));
%        uuu=[uuu derexactau-Sderu];
%       vvv=[vvv derexactav-Sderv];
%       vvv=[vvv exacta-zz];
      error2=(exacta-zz).^2;
      %dererror2=(derexactau-Sderu).^2 + (derexactav-Sderv).^2;

%
% %       pause
%       vvv=[vvv derexactav-Sderv];
%       pause
%       dererror2
%       pause

      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      L2_error = L2_error + sum(sum(error2)) * J2 * wt;
     % L2_dererror = L2_dererror + sum(sum(dererror2)) * J2 * wt;
%       uuu=[uuu sum(sum(error2)) * J2 * wt];
    end
end
%H1_error=sqrt(L2_error + L2_dererror)
L2_error=sqrt(L2_error)
