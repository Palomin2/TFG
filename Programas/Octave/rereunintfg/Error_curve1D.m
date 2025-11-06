% U=u;
L2_error=0;
L2_dererror=0;
clear exacta
k=1;
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
      clear zz
      clear derexacta
      clear aux2
      clear dxx
      clear dyy

      xx=zeros(1,length(Xi));
      yy=zeros(1,length(Xi));
      zz=zeros(1,length(Xi));
      dxx=zeros(1,length(Xi));
      dyy=zeros(1,length(Xi));

      exacta=zeros(1,length(Xi));
      derexacta=zeros(1,length(Xi));

      s_gp = interp1(xi_samp, s_samp, Xi, 'cubic');
      % s_gp=Xi;
      for j=1:length(Xi)
        C=NURBSCurvepoint3D(length(U),p,knotVec,[controlPts(:,1),controlPts(:,2),U,controlPts(:,3)],Xi(j));

        % NO ME FUNCIONA EL ERROR EN H1, LA DERIVADA SE ESCAPA...

        D=dersNURBSCurvepoint3D(length(U),p,knotVec,[controlPts(:,1),controlPts(:,2),U,controlPts(:,3)],Xi(j)); % Hacer que funcione esto?
        dzz(j)=D(3);

        % D=dersNURBSCurvepoint(length(U),p,knotVec,[controlPts(:,1),U,controlPts(:,3)],Xi(j)); % Como esto funciona en el caso de una recta, puede ser que no haya convergencia de la derivada desde el punto de vista teorico
        % dzz(j)=D(2);

        % if Xi>=0 & Xi < 1
        %     D1=NURBSCurvepoint3D(length(U),p,knotVec,[controlPts(:,1),controlPts(:,2),U,controlPts(:,3)],Xi(j)+1e-6); %
        %     D2=NURBSCurvepoint3D(length(U),p,knotVec,[controlPts(:,1),controlPts(:,2),U,controlPts(:,3)],Xi(j)); %
        %     D=(D1-D2)/1e-6;
        % elseif  Xi == 1
        %     D1=NURBSCurvepoint3D(length(U),p,knotVec,[controlPts(:,1),controlPts(:,2),U,controlPts(:,3)],Xi(j)); %
        %     D2=NURBSCurvepoint3D(length(U),p,knotVec,[controlPts(:,1),controlPts(:,2),U,controlPts(:,3)],Xi(j)-1e-6); %
        %     D=(D1-D2)/1e-6;
        % end
        % dzz(j)=D(3);

        xx(j)=C(1);
        yy(j)=C(2);
        zz(j)=C(3);

        exacta(j)=sin(c*pi*s_gp(j));
        derexacta(j)=c*pi*cos(c*pi*s_gp(j));

        Z(k)=dzz(j)-derexacta(j);
        k=k+1;
      end
      error2=(exacta-zz).^2;
      dererror2=(derexacta-dzz).^2;

      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      L2_error = L2_error + sum(error2) * J2 * wt;
      L2_dererror = L2_dererror + sum(dererror2) * J2 * wt;
    end
end
H1_error=sqrt(L2_error + L2_dererror) % HAY CONVERGENCIA EN L2?
L2_error=sqrt(L2_error)


% Probar un ejemplo con B-splines?
% plot(Z)
