% U=u;
L2_error=0;
L2_dererror=0;
clear exacta
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
      clear exacta
      clear aux
      clear xx
      clear yy
      clear derexacta
      clear aux2
      clear dxx
      clear dyy
      xx=zeros(1,length(Xi));
      yy=zeros(1,length(Xi));
      exacta=zeros(1,length(Xi));
      dxx=zeros(1,length(Xi));
      dyy=zeros(1,length(Xi));
      derexacta=zeros(1,length(Xi));
      for j=1:length(Xi)
        C=NURBSCurvepoint(length(U),p,knotVec,[controlPts(:,1),U,controlPts(:,3)],Xi(j));
        D=dersNURBSCurvepoint(length(U),p,knotVec,[controlPts(:,1),U,controlPts(:,3)],Xi(j));
%         D
%         pause
        aux(j,:)=C;
        aux2(j,:)=D;
        xx(j)=C(1);
        yy(j)=C(2);
        dxx(j)=D(1);
        dyy(j)=D(2);
        exacta(j)=sin(c*pi*xx(j));
        derexacta(j)=c*pi*cos(c*pi*xx(j));
%         derexacta(j)-c*pi*cos(c*pi*Xi(j))
%         pause
%         dxx(j)
%         derexacta(j)
%         dyy(j)
%         derexacta(j)-dyy(j)
%         pause
      end
      error2=(exacta-yy).^2;
      
      dererror2=(derexacta-dyy/dxx).^2;
      
      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      L2_error = L2_error + sum(error2) * J2 * wt;
      L2_dererror = L2_dererror + sum(dererror2) * J2 * wt;
    end
end
H1_error=sqrt(L2_error + L2_dererror)
L2_error=sqrt(L2_error)