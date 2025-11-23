% Vinh Phu Nguyen
% Delft University of Technology
% Adapted from ISOGAT, A.V. Vuong.

for c=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
% %     uKnotVectorUB = unique(uBKnot);
% %     uKnotVectorVB = unique(vBKnot);
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
% %     newKnotsXB = uKnotVectorUB(1:end-1) + 0.5*diff(uKnotVectorUB);
% %     newKnotsYB = uKnotVectorVB(1:end-1) + 0.5*diff(uKnotVectorVB);
    %% h-refinement (NURBS) in x-direction
    dim = size(controlPts,2); % numero de columnas de la matriz controlPts
    
    nonewkX      = size(newKnotsX,2);
    newprojcoord = zeros(noPtsX*noPtsY+nonewkX*noPtsY,dim+1);
    
% %     nonewkXB      = size(newKnotsXB,2);
% %     newprojcoordB = zeros(nB*mB+nonewkXB*mB,3);
    
    rstart = 1;
    wstart = 1;
    
    for j=1:noPtsY
        rstop = rstart + noPtsX-1;
        wstop = wstart + noPtsX-1 + nonewkX;
        
        locCP        = controlPts(rstart:rstop,:);
        locweights   = weights   (rstart:rstop);
        locprojcoord = nurb2proj(noPtsX, locCP, locweights);
        
        % refinement of x
        [tempknotVectorX,tempControlPts] = ...
            RefineKnotVectCurve(noPtsX-1,p,uKnot,locprojcoord,newKnotsX,nonewkX-1);
        
        newprojcoord(wstart:wstop,:)=tempControlPts;
        wstart = wstop+1;
        rstart = rstop+1;
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MIO:
% %         for k=1:mB
% %             rcpindexB     = k:mB:nB*mB;
% %             locCPB        = Bij(rcpindexB,:);
% %             locweightsB   = weightsB(rcpindexB);
% %             locprojcoordB = nurb2proj(nB, locCPB, locweightsB);
% %             
% %             [tempknotVectorBkX, Buk]= ...
% %                 RefineKnotVectCurve(nB-1,pB,uBKnot,locprojcoordB,newKnotsXB,nonewkXB-1);
% %             
% % %             newB(:,k,1)=Buk(:,1);
% % %             newB(:,k,2)=Buk(:,2);
% %             wcpindexB = k:mB:mB*(nB+nonewkXB);
% %             newprojcoordB(wcpindexB,:)=Buk;
% %         end
% %         clear B
% %         %clear PBw
% %         clear Buk
%         B(:,:,:)=newB(:,:,:);
%         clear newB
        
% %         uBKnot = tempknotVectorBkX;
% %         [Bij, weightsB]=proj2nurbs(newprojcoordB);
% %         nB = nB+nonewkXB;
        
        
    uKnot                 = tempknotVectorX;
    [controlPts, weights] = proj2nurbs(newprojcoord);
    noPtsX                = noPtsX+nonewkX;
    %% h-refinement (NURBS) in y-direction)
    nonewkY       = size(newKnotsY,2);
    newprojcoord  = zeros(noPtsX*noPtsY+nonewkY*noPtsX,dim+1);
    
% %     nonewkYB       = size(newKnotsYB,2);
% %     newprojcoordB  = zeros(nB*mB+nonewkYB*nB,3);
       
    for i=1:noPtsX
        %create index for reading controlPoints
        rcpindex     = i:noPtsX:noPtsX*noPtsY;
        locCP        = controlPts(rcpindex,:);
        locweights   = weights(rcpindex);
        locprojcoord = nurb2proj(noPtsY, locCP, locweights);
        
        % refinement of y
        [tempknotVectorY,tempcontrolPoints] = ...
            RefineKnotVectCurve(noPtsY-1,q,vKnot,locprojcoord,newKnotsY,nonewkY-1);
        
        wcpindex = i:noPtsX:noPtsX*(noPtsY+nonewkY);
        newprojcoord(wcpindex,:) = tempcontrolPoints;
    end
      %%%%%%%%%%%%%%%%%%%%%MIO:
% %       rstartB = 1;
% %       wstartB = 1;
% %         for k=1:nB
% %         rstopB = rstartB + mB-1;
% %         wstopB = wstartB + mB-1 + nonewkYB;
% %         
% %         locCPB        = Bij(rstartB:rstopB,:);
% %         locweightsB   = weightsB(rstartB:rstopB);
% %         locprojcoordB = nurb2proj(mB, locCPB, locweightsB);
% %             [tempknotVectorBkY, Bvk]= ...
% %                 RefineKnotVectCurve(mB-1,qB,vBKnot,locprojcoordB,newKnotsYB,nonewkYB-1);
% % %             newB(k,:,1)=Bvk(:,1);
% % %             newB(k,:,2)=Bvk(:,2);
% %         newprojcoordB(wstartB:wstopB,:)=Bvk;
% %         wstartB = wstopB+1;
% %         rstartB = rstopB+1;
% %         end
% % %         clear B
% %         clear Bvk
% % %         B(:,:,:)=newB(:,:,:);
% % %         clear newB
% %         vBKnot = tempknotVectorBkY;
% %         [Bij, weightsB] = proj2nurbs(newprojcoordB);
% %         mB = mB + nonewkYB;
        
    vKnot                 = tempknotVectorY;
    [controlPts, weights] = proj2nurbs(newprojcoord);
    noPtsY                = noPtsY + nonewkY;
    WI{c+1}=weights;
end