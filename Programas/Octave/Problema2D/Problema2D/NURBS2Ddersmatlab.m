function [dRdxi, dRdeta] = NURBS2Ddersmatlab(Xi,Eta,p,q,uKnot,vKnot,weights)
n=length(uKnot)-p-1;
m=length(vKnot)-q-1;
w=0;
dwdxi=0;
dwdeta=0;
dRdxi=zeros(1,(p+1)*(q+1));
dRdeta=zeros(1,(p+1)*(q+1));
uspan=findspan(n,p,Xi,uKnot); 
vspan=findspan(m,q,Eta,vKnot);
N = BasisFuns(uspan+1,Xi,p,uKnot);
M = BasisFuns(vspan+1,Eta,q,vKnot);
uders = DerBasisFuns(uspan+1,Xi,p,n,uKnot); 
vders = DerBasisFuns(vspan+1,Eta,q,m,vKnot); 

uind=uspan-p;
for j=1:q+1
    vind=vspan-q+j-1;
    for i=1:p+1
        c=uind+i+n*vind;
        wgt=weights(c);
        w=w+N(i)*M(j)*wgt;
        dwdxi=dwdxi+uders(2,i)*M(j)*wgt;
        dwdeta=dwdeta+vders(2,j)*N(i)*wgt;
    end
end 
k=1;
for j=1:q+1
    vind=vspan-q+j-1;
    for i=1:p+1
        c=uind+i+n*vind;
        fac= weights(c)/(w*w);
        % R(k)=N(i)*M(j)*fac*w;
        dRdxi(k)=(uders(2,i)*M(j)*w - N(i)*M(j)*dwdxi)*fac;
        dRdeta(k)=(vders(2,j)*N(i)*w - N(i)*M(j)*dwdeta)*fac;
        k=k+1;
    end
end 
