function C = dersNURBSCurvepoint3D(n,p,U,Pw,u0)
span=findspan(n-1,p,u0,U);
N = BasisFuns(span+1,u0,p,U); 
D = DerBasisFuns(span+1,u0,p,n,U); 

w=0;
w=sum(N.*Pw(span+1-p:span+1,4)');
derw= sum(D(2,:).*Pw(span+1-p:span+1,4)'); 
dRdxi=((D(2,:).*Pw(span+1-p:span+1,4)')*w - derw*(N.*Pw(span+1-p:span+1,4)'))/w^2;

dRdx=dRdxi;

Cw=zeros(1,4);
for i=1:p+1
    Cw=Cw+dRdx(i)*Pw(span-p+i,:);
%     Pw(span-p+i,:)
%     Cw
%     pause
%     Cw=Cw+D(2,i)*Pw(span-p+i,:);
end 
% Cw
% pause
% Cw(1,[1,2,3])=Cw(1,[1,2,3])./Cw(1,[4,4,4]);

C=Cw;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % function dC = dersNURBSCurvepoint3D(n,p,U,Pw,u0)
% % % % % Devuelve la derivada cartesiana dC/dxi de la curva NURBS 3D en u0.
% % % % % n = number of control points
% % % % % Pw = n x 4 matrix: [x y z w]
% % % % % Resultado: dC (1x3)
% % % % 
% % % % % seguridad en extremos: evitar evaluar exactamente en U(end)
% % % % tol = 1e-12;
% % % % if abs(u0 - U(end)) < tol
% % % %     u0 = U(end) - tol;
% % % % elseif abs(u0 - U(1)) < tol
% % % %     u0 = U(1) + tol;
% % % % end
% % % % 
% % % % % span y bases
% % % % span = findspan(n-1, p, u0, U);           % convención tuya: n-1 como last index
% % % % idx  = (span+1-p) : (span+1);             % índices locales para N, D y Pw
% % % % 
% % % % N   = BasisFuns(span+1, u0, p, U);        % 1 x (p+1)
% % % % D1  = DerBasisFuns(span+1, u0, p, 1, U);  % (nd+1) x (p+1) -> D1(2,:) es N'
% % % % 
% % % % % obtener pesos locales
% % % % w_i = Pw(idx,4)';                         % 1 x (p+1)
% % % % 
% % % % % W y W' (escalares)
% % % % W  = sum(N .* w_i);
% % % % Wd = sum(D1(2,:) .* w_i);
% % % % 
% % % % % derivada racional R'_i(ξ)
% % % % % cuidado con formas: (1 x (p+1)) .* (1 x (p+1)) => 1 x (p+1)
% % % % dRdxi = ( (D1(2,:) .* w_i) * W - Wd * (N .* w_i) ) / (W^2);  % 1 x (p+1)
% % % % 
% % % % % derivada cartesiana: dC/dxi = sum_i dRdxi(i) * P_i (P_i = Pw(idx,1:3) are Cartesian)
% % % % Pcoords = Pw(idx,1:3);                    % (p+1) x 3
% % % % dC = dRdxi * Pcoords;                     % 1 x 3  (row vector)
% % % % 
% % % % end