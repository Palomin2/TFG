function [Ap , stern]=Api(N,p,i,h)
q=4; % si pones global q, tarda mucho!
Ap=0;
if q==2
    stern=stencil2(i,N,h);
    if (i>2 && i<N-1)
        Ap=stern*p(i-q:i+q);
    elseif i==1
        Ap=stern(3)*p(i);
    elseif i==2
        Ap=stern(2:5)*p(i-1:i+2);
    elseif i==N-1
        Ap=stern(1:4)*p(i-2:i+1);
    elseif i==N
        Ap=stern(3)*p(i);
    end
    
elseif q==3
    stern=stencil3(i,N,h);
    if i==1
        Ap=stern(4)*p(i);
    elseif i==2
        Ap=stern(3:7)*p(i-1:i+3);
    elseif i==3
        Ap=stern(2:7)*p(i-2:i+3);
    elseif i==N-2
        Ap=stern(1:6)*p(i-3:i+2);
    elseif i==N-1
        Ap=stern(1:5)*p(i-3:i+1);
    elseif i==N
        Ap=stern(4)*p(i);
    elseif (i>3 && i<N-2)
        Ap=stern*p(i-q:i+q);
    end
elseif q==4
    stern=stencil4(i,N,h);
    if i==1
        Ap=stern(5)*p(i);
    elseif i==2
        Ap=stern(4:9)*p(i-1:i+4);
    elseif i==3
        Ap=stern(3:9)*p(i-2:i+4);
    elseif i==4
        Ap=stern(2:9)*p(i-3:i+4);
    elseif i==N-3
        Ap=stern(1:8)*p(i-4:i+3);
    elseif i==N-2
        Ap=stern(1:7)*p(i-4:i+2);
    elseif i==N-1
        Ap=stern(1:6)*p(i-4:i+1);
    elseif i==N
        Ap=stern(5)*p(i);
    elseif (i>4 && i<N-3)
        Ap=stern*p(i-q:i+q);
    end
elseif q==5
    stern=stencil5(i,N,h);
    if i==1
        Ap=stern(6)*p(i);
    elseif i==2
        Ap=stern(5:11)*p(i-1:i+5);
    elseif i==3
        Ap=stern(4:11)*p(i-2:i+5);
    elseif i==4
        Ap=stern(3:11)*p(i-3:i+5);
    elseif i==5
        Ap=stern(2:11)*p(i-4:i+5);
    elseif i==N-4
        Ap=stern(1:10)*p(i-5:i+4);
    elseif i==N-3
        Ap=stern(1:9)*p(i-5:i+3);
    elseif i==N-2
        Ap=stern(1:8)*p(i-5:i+2);
    elseif i==N-1
        Ap=stern(1:7)*p(i-5:i+1);
    elseif i==N
        Ap=stern(6)*p(i);
    elseif (i>5 & i<N-4)
        Ap=stern*p(i-5:i+5);
    end
elseif q==6
    stern=stencil6(i,N,h);
    if (i>6 & i<N-5)
        Ap=stern*p(i-6:i+6);
    elseif i==1
        Ap=stern(7)*p(i);
    elseif i==2
        Ap=stern(6:13)*p(i-1:i+6);
    elseif i==3
        Ap=stern(5:13)*p(i-2:i+6);
    elseif i==4
        Ap=stern(4:13)*p(i-3:i+6);
    elseif i==5
        Ap=stern(3:13)*p(i-4:i+6);
    elseif i==6
        Ap=stern(2:13)*p(i-5:i+6);
    elseif i==N-5
        Ap=stern(1:12)*p(i-6:i+5);
    elseif i==N-4
        Ap=stern(1:11)*p(i-6:i+4);
    elseif i==N-3
        Ap=stern(1:10)*p(i-6:i+3);
    elseif i==N-2
        Ap=stern(1:9)*p(i-6:i+2);
    elseif i==N-1
        Ap=stern(1:8)*p(i-6:i+1);
    elseif i==N
        Ap=stern(7)*p(i);
    end
elseif q==7
    stern=stencil7(i,N,h);
    if i==1
        Ap=stern(8)*p(i);
    elseif i==2
        Ap=stern(7:15)*p(i-1:i+7);
    elseif i==3
        Ap=stern(6:15)*p(i-2:i+7);
    elseif i==4
        Ap=stern(5:15)*p(i-3:i+7);
    elseif i==5
        Ap=stern(4:15)*p(i-4:i+7);
    elseif i==6
        Ap=stern(3:15)*p(i-5:i+7);
    elseif i==7
        Ap=stern(2:15)*p(i-6:i+7);
    elseif i==N-6
        Ap=stern(1:14)*p(i-7:i+6);
    elseif i==N-5
        Ap=stern(1:13)*p(i-7:i+5);
    elseif i==N-4
        Ap=stern(1:12)*p(i-7:i+4);
    elseif i==N-3
        Ap=stern(1:11)*p(i-7:i+3);
    elseif i==N-2
        Ap=stern(1:10)*p(i-7:i+2);
    elseif i==N-1
        Ap=stern(1:9)*p(i-7:i+1);
    elseif i==N
        Ap=stern(8)*p(i);
    elseif (i>7 & i<N-6)
        Ap=stern*p(i-7:i+7);
    end
elseif q==8
    stern=stencil8(i,N,h);
    if i==1
        Ap=stern(9)*p(i);
    elseif i==2
        Ap=stern(8:17)*p(i-1:i+8);
    elseif i==3
        Ap=stern(7:17)*p(i-2:i+8);
    elseif i==4
        Ap=stern(6:17)*p(i-3:i+8);
    elseif i==5
        Ap=stern(5:17)*p(i-4:i+8);
    elseif i==6
        Ap=stern(4:17)*p(i-5:i+8);
    elseif i==7
        Ap=stern(3:17)*p(i-6:i+8);
    elseif i==8
        Ap=stern(2:17)*p(i-7:i+8);
    elseif i==N-7
        Ap=stern(1:16)*p(i-8:i+7);
    elseif i==N-6
        Ap=stern(1:15)*p(i-8:i+6);
    elseif i==N-5
        Ap=stern(1:14)*p(i-8:i+5);
    elseif i==N-4
        Ap=stern(1:13)*p(i-8:i+4);
    elseif i==N-3
        Ap=stern(1:12)*p(i-8:i+3);
    elseif i==N-2
        Ap=stern(1:11)*p(i-8:i+2);
    elseif i==N-1
        Ap=stern(1:10)*p(i-8:i+1);
    elseif i==N
        Ap=stern(9)*p(i);
    elseif (i>8 & i<N-7)
        Ap=stern*p(i-8:i+8);
    end
end
