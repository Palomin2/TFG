function [pf,n,uKnotf,qf,m,vKnotf,Qf,Wf] = prefinement(p,noPtsX,uKnot,q,noPtsY,vKnot,P,weights)
    
pf=p+1;
qf=q+1;
n=noPtsX+1;
m=noPtsY+1;
uKnotf=[0 uKnot 1];
vKnotf=[0 vKnot 1];
matriz=zeros(n*m,n*m);
for ii=1:n
    y(ii)=(ii-1)/(n-1); 
    for jj=1:m
        x(jj)=(jj-1)/(m-1); 
        ind1=m*(ii-1)+jj;
        for i=1:n
            for j=1:m
                ind2=m*(i-1)+j;
                matriz(ind1,ind2)=OneBasisFun(p+1,length(uKnotf),uKnotf,j,x(jj))*OneBasisFun(p+1,length(vKnotf),vKnotf,i,y(ii)) ;
            end
        end
    end
end
    
matrizrhs=zeros(n*m,noPtsX*noPtsY);
for ii=1:n
    y(ii)=(ii-1)/(n-1); 
    for jj=1:m
        x(jj)=(jj-1)/(m-1); 
        ind1=m*(ii-1)+jj;
        for i=1:noPtsY
            for j=1:noPtsX
                ind2=noPtsX*(i-1)+j;
                matrizrhs(ind1,ind2)=OneBasisFun(p,length(uKnot),uKnot,j,x(jj))*OneBasisFun(p,length(vKnot),vKnot,i,y(ii)) ;
            end
        end
    end
end

rhsx=matrizrhs*(P(:,1).*weights);
rhsy=matrizrhs*(P(:,2).*weights);
rhsw=matrizrhs*weights;
X=matriz\rhsx;
Y=matriz\rhsy;
Wf=matriz\rhsw;
Qf=zeros(n*m,2);
Qf(:,1)=X./Wf;
Qf(:,2)=Y./Wf;