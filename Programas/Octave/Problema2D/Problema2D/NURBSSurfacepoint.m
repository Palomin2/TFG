function S = NURBSSurfacepoint(n,p,U,m,q,V,P,u0,v0)
uspan=findspan(n-1,p,u0,U); 
vspan=findspan(m-1,q,v0,V);
Pw=zeros(size(P));
Sw=zeros(4,1);
temp=zeros(4,q+1);
Nu = BasisFuns(uspan+1,u0,p,U);
Nv = BasisFuns(vspan+1,v0,q,V);

Pw=P;
Pw(:,:,[1,2,3])=Pw(:,:,[1,2,3]).*Pw(:,:,[4,4,4]);
aux=zeros(4,1);
for i=1:q+1
    for j=1:p+1
        %aux=zeros(4,1);
        for l=1:4
            aux(l,1)=Nu(j)*Pw(uspan-p+j,vspan-q+i,l);
        end
        temp(:,i)=temp(:,i)+aux;
        %temp(:,i)=temp(:,i)+Nu(j)*Pw(uspan-p+j,vspan-q+i,:);
    end
end
for k=1:q+1
    Sw(:)=Sw(:)+Nv(k)*temp(:,k);
end
% Sw
% pause
Sw([1,2,3],1)=Sw([1,2,3],1)./Sw([4,4,4],1);
% Sw
% pause
S=Sw;
% 
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:51
%   for j=1:51 
%       %S=NURBSSurfacepoint(noPtsX,p,uKnot,noPtsY,q,vKnot,P,2*(i-1)/50,2*(j-1)/50);
%       S=NURBSSurfacepoint(noPtsX,p,uKnot,noPtsY,q,vKnot,P,1.8,1.8);
%       xx(i)=S(1);
%       yy(j)=S(2);
%       zz(i,j)=S(3);
%       %ww(i,j)= xx(i)^2*yy(j)^2 -2*(xx(i)^2*yy(j)) -2*(xx(i)*yy(j)^2) +4*xx(i)*yy(j);
%       ww(i,j)=sin(pi*(S(1))/2)*sin(pi*(S(2))/2);
%   end
% end 



%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


% 
% uspan=findspan(noPtsX-1,p,1.8,uKnot); 
% vspan=findspan(noPtsY-1,q,1.8,vKnot);
% Pw=zeros(size(P));
% Sw=zeros(4,1);
% temp=zeros(4,q+1);
% Nu = BasisFuns(uspan+1,1.8,p,uKnot); 
% Nv = BasisFuns(vspan+1,1.8,q,vKnot); 
% 
% Pw=P;
% Pw(:,:,[1,2,3])=Pw(:,:,[1,2,3]).*Pw(:,:,[4,4,4]);
% 
% aux=zeros(4,1);
% for i=1:q+1
%     for j=1:p+1
%         for l=1:4
%             aux(l,1)=Nu(j)*Pw(uspan-p+j,vspan-q+i,l);
%         end
%         temp(:,i)=temp(:,i)+aux;
%         %temp(:,i)=temp(:,i)+Nu(j)*Pw(uspan-p+j,vspan-q+i,:);
%     end
% end
% for k=1:q+1
%     Sw(:)=Sw(:)+Nv(k)*temp(:,k);
% end
% Sw([1,2,3],1)=Sw([1,2,3],1)./Sw([4,4,4],1);
% %Sw(4)
% S=Sw;