gamma=zeros(1,p+2);
for i=1:p+2
    gamma(i)=(i-1)/(p+1);
end
a = 0.3; % inner radius
b = 0.5; % outer radius
controlPtsref=[];
for j=1:noPtsY
    P=controlPts((j-1)*noPtsX+1:j*noPtsX,:);
    %P=[a 0; 0.5*(a+b) 0 ; b 0];
    Pweights=weights((j-1)*noPtsX+1:j*noPtsX);
    CP=[P, Pweights];
    newCP=zeros(p+2,3);
    newCP(1,:)=CP(1,:);
    newCP(noPtsX+1,:)=CP(noPtsX,:);
    for i=2:noPtsX
        newCP(i,:)=CP(i,:)*(1-gamma(i))+CP(i-1,:)*gamma(i);
    end
    %[A1,A2]=bspdegelev(2,[P';weights'],[0 0 0 1 1 1],1);
    newCP
    pause
    controlPtsref=[controlPtsref; newCP]
end 
    
    
    
    