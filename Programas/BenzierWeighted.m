%Benzier weighted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constantes de control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=100; %numero de pasos
step=1/m;
u=0:step:1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINIMOS LOS P (puntos de control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for w=0:0.1:20
  P=[0, 1, 2, 3, 1.5, 0;
     0, 5, 5, 0, -3, 0];
  we=[1;1*w*3;1+w*2;1+7*w; 1+2*w; 1];
  % we=we/norm(w,1);

  [aux n] = size(P); %guardamos el grado n

  sol=zeros(aux,m+1);

  for k=1:m+1
  coef=0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %DEFINIMOS B(i,n) por recursividad
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Bold=zeros(n);
    Bnew=zeros(n);

    Bold(1)=1; %Caso B(0,0) n=0

    for j=2:n  %bucle para el indice i
      Bnew=zeros(n);
      for i=1:j-1%bucle para el indice n
          Bnew(i)=Bnew(i) + (1-u(k))*Bold(i);
          Bnew(i+1)=Bnew(i+1) + (u(k))*Bold(i);
      endfor
      Bold=Bnew;
    endfor
    for j=1:n
      coef= coef+Bnew(j)*we(j);
    endfor
    for j=1:n
      sol(:,k)=sol(:,k)+(Bnew(j)*P(:,j)*we(j))/coef;
    endfor
  endfor

  plot(sol(1,:),sol(2,:), color='b',P(1,:),P(2,:), color = 'r');
  xlim([0 3])
  ylim([-5 5])
  grid on
  pause(0.01)
endfor


