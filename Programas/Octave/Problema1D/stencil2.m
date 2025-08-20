function stern = stencil2(i,N,h)
if i==1 %longitud 5
    stern=[0 0 4/3 0 0];
elseif i==2
    stern=[0,0,4/3,-1/6,-1/6];
elseif i==3
    stern=[0,-1/6,1,-1/3,-1/6];
elseif i==N-2
    stern=[-1/6,-1/3,1,-1/6,0];
elseif i==N-1
    stern=[-1/6,-1/6,4/3,0,0];
elseif i==N
    stern=[0 0 4/3 0 0];
else
    stern=[-1/6,-1/3,1,-1/3,-1/6];
end
                                                                                                                                                                                                                                                                                                                            
stern=stern/h;