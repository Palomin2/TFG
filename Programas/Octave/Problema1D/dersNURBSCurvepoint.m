function C = dersNURBSCurvepoint(n,p,U,Pw,u0)
span=findspan(n-1,p,u0,U);
N = BasisFuns(span+1,u0,p,U); 
D = DerBasisFuns(span+1,u0,p,n,U); 
% Pw(span+1-p:span+1,:)
% pause
w=0;
w=sum(BasisFuns(span+1,u0,p,U).*Pw(span+1-p:span+1,3)');
derw= sum(D(2,:).*Pw(span+1-p:span+1,3)'); 
dRdxi=((D(2,:).*Pw(span+1-p:span+1,3)')*w - derw*(N.*Pw(span+1-p:span+1,3)'))/w^2;
% jacob1 = dRdxi*Pw(span+1-p:span+1,1:2);
% J1     = norm (jacob1);
% dRdx   = (1/J1)*dRdxi;
dRdx=dRdxi;
% Pw(:,[1,2])=Pw(:,[1,2]).*Pw(:,[3,3]);
Cw=zeros(1,3);
for i=1:p+1
    Cw=Cw+dRdx(i)*Pw(span-p+i,:);
%     Pw(span-p+i,:)
%     Cw
%     pause
%     Cw=Cw+D(2,i)*Pw(span-p+i,:);
end 
% Cw
% pause
% Cw(1,[1,2])=Cw(1,[1,2])./Cw(1,[3,3]);
C=Cw;
% n=5
% p=2
% Pw = [0,0,1 ; 1,1,1; 3,2,1; 4,1,1; 5,-1,1]
% U = [0,0,0,1,2,3,3,3]
% 
% for j=1:200
%     C=NURBSCurvepoint(n,p,U,Pw,3*(j-1)/200);
%     x(j)=C(1); 
%     y(j)=C(2); 
% end 
% x(201)=5;y(201)=-1;
% plot(x,y,Pw(:,1),Pw(:,2),'black');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% aux=zeros(200,3);
% for j=1:201
%     C=NURBSCurvepoint(length(U),p,knotVec,[controlPts(:,1),U,controlPts(:,3)],(j-1)/200);
%     aux(j,:)=C;
%     xx(j)=C(1);
%     yy(j)=C(2); % solucion
% end 
%  %xx(201)=1;yy(201)=0;
% plot(xx,yy);
% 
% abcisas = 0:1/200:1 ;
% ordenadas = (-abcisas.^3 + abcisas)/6 ;
% ordenadas = sin(pi*xx) ;
% plot(controlPts(:,1),U,'black',xx,ordenadas,'red',xx,yy,'blue')
% 
% max(abs(yy-ordenadas))
% norm(yy-ordenadas,2)
% 
% % GRAFICA DEL ERROR
% 
% plot(xx,yy-ordenadas)
