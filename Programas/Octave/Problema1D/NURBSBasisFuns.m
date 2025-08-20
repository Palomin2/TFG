function N = NURBSBasisFuns(number,i,x,p,U,w) 
aux=BasisFuns(i,x,p,U).*w([i-p:i]);
N=aux(number)/sum(aux);

% %% Pintar una funcion base:
% n=7
% w=[1,1,1,1,1,1,1]
% p=3
% U = [0,0,0,0,1/4,1/2,3/4,1,1,1,1]
% %% i=findspan(n,p,x,U)+1;
% 
% for j=1:201
%     x(j)=(j-1)/200; %% particion (0,2)
%     if findspan(n,p,x(j),U)==p
%       y1(j)=NURBSBasisFuns(1,p+1,x(j),p,U,w); % solucion
%     end
% end 
% plot(x(1:length(y1)),y1);
% 
% %%%%%%%%%%%%%%%%% PINTAR TODAS LAS FUNCIONES BASE
% for j=1:201
%     x(j)=(j-1)/200; %% particion (0,2)
%     if findspan(n,p,x(j),U)==p
%       y1(j)=NURBSBasisFuns(1,p+1,x(j),p,U,w); % solucion
%       y2(j)=NURBSBasisFuns(2,p+1,x(j),p,U,w); % solucion
%       y3(j)=NURBSBasisFuns(3,p+1,x(j),p,U,w); % solucion
%       y4(j)=NURBSBasisFuns(4,p+1,x(j),p,U,w); % solucion
%     end
%     if findspan(n,p,x(j),U)==p+1
%       y2(j)=NURBSBasisFuns(1,p+2,x(j),p,U,w); % solucion
%       y3(j)=NURBSBasisFuns(2,p+2,x(j),p,U,w); % solucion
%       y4(j)=NURBSBasisFuns(3,p+2,x(j),p,U,w); % solucion
%       y5(j-50)=NURBSBasisFuns(4,p+2,x(j),p,U,w); % solucion
%     end
%     if findspan(n,p,x(j),U)==p+2
%       y3(j)=NURBSBasisFuns(1,p+3,x(j),p,U,w); % solucion
%       y4(j)=NURBSBasisFuns(2,p+3,x(j),p,U,w); % solucion
%       y5(j-50)=NURBSBasisFuns(3,p+3,x(j),p,U,w); % solucion
%       y6(j-100)=NURBSBasisFuns(4,p+3,x(j),p,U,w); % solucion
%     end
%     if findspan(n,p,x(j),U)==p+3
%       y4(j)=NURBSBasisFuns(1,p+4,x(j),p,U,w); % solucion
%       y5(j-50)=NURBSBasisFuns(2,p+4,x(j),p,U,w); % solucion
%       y6(j-100)=NURBSBasisFuns(3,p+4,x(j),p,U,w); % solucion
%       y7(j-150)=NURBSBasisFuns(4,p+4,x(j),p,U,w); % solucion
%     end
%     if findspan(n,p,x(j),U)==p+4
%       y5(j-50)=NURBSBasisFuns(1,p+5,x(j),p,U,w); % solucion
%       y6(j-100)=NURBSBasisFuns(2,p+5,x(j),p,U,w); % solucion
%       y7(j-150)=NURBSBasisFuns(3,p+5,x(j),p,U,w); % solucion
%     end
%     if findspan(n,p,x(j),U)==p+5
%       y6(j-100)=NURBSBasisFuns(1,p+6,x(j),p,U,w); % solucion
%       y7(j-150)=NURBSBasisFuns(2,p+6,x(j),p,U,w); % solucion
%     end
%     if findspan(n,p,x(j),U)==p+6
%       y7(j-150)=NURBSBasisFuns(1,p+7,x(j),p,U,w); % solucion
%     end
% end 
% y5(151)=0;
% y6(101)=0;
% y7(51)=1;
% plot(x(1:length(y1)),y1,x(1:length(y2)),y2,x(1:length(y3)),y3,x(1:length(y4)),y4,x(length(y1)+1:201),y5,x(length(y2)+1:201),y6,x(length(y3)+1:201),y7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
