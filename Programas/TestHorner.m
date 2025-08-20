m=10; %numero de pasos
a=0; %valor inicial del param de la curva
b=2; %valor final del param de la curva
step=1/m; %valor del pasos
u=a:step:b; %vector del param discreto con paso step


for k=0:100
  a=[0,1,1;  %primera coordenada (x)
     0,-k/10,k/10];  %segunda coordenada (y)
     %vector de puntos
  [aux n] = size(a);



  %C=zeros(aux,n);
  curva=zeros(aux,m);



  for i=1:m

    C = a(:,n);

    for j=n-1:-1:1
      C = C*u(i) + a(:,j);

    endfor
    curva(:,i)=C;

  endfor

  x=zeros(aux,n);

  x(:,1)=a(:,1);
  for j=2:n
    j;
    x(:,j)=a(:,j)+x(:,j-1);
  endfor
  %hold on
  %no funciona?

  plot(curva(1,:),curva(2,:), color='b',a(1,:),a(2,:), color = 'r*',
  x(1,:),x(2,:), color='g*');
  xlim([0 5])
  ylim([-5 5])
  grid on
  %plot(a(1,:),a(2,:), color = 'r*');
  %xlim=([0 3])
  %ylim=([0 3])
  %hold off
  pause(0.1)
endfor
