%Benzier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constantes de control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=100; %numero de pasos
u=0:1/m:1; %vector del parametro de la curva

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINIMOS LOS P (puntos de control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


P=[0, 1, 2, 3;
     0, 5, 0, 4];

[aux n] = size(P); %guardamos el grado n (tener en cuenta que sera n+1 realmente
                   % ya que octave indexiza desde el 1)

sol=zeros(aux,m+1); %sera nuestra curva
sol2=zeros(aux,m+1);%sera el vector la recta tangente (solo para n=3)

for k=1:m+1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %DEFINIMOS B(i,n) por recursividad
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Bold=zeros(n);%reinicio de las variables para cada punto de la curva
  %Bnew=zeros(n); &no es necesaria pero la deja comentada por claridad

  Bold(1)=1; %Caso B(0,0) n=0

  for j=2:n
    Bnew=zeros(n);
    for i=1:j-1
      Bnew(i)=Bnew(i) + (1-u(k))*Bold(i);
      Bnew(i+1)=Bnew(i+1) + (u(k))*Bold(i);
    endfor
    Bold=Bnew;
  endfor

  %declaracion de la curva
  for j=1:n
    sol(:,k)=sol(:,k)+Bnew(j)*P(:,j);
  endfor

endfor

for uo=0:0.01:1 %bucle que mueve la tangente y el punto tg
    for k=1:m+1
     sol2(:,k)=(P(:,2)*uo+(1-uo)*P(:,1))*(1-u(k))+(P(:,3)*uo+(1-uo)*P(:,2))*(u(k));
    endfor

    ptoTg=(P(:,2)*uo+(1-uo)*P(:,1))*(1-uo)+(P(:,3)*uo+(1-uo)*P(:,2))*(uo);

  plot(sol(1,:),sol(2,:), color='b',P(1,:),P(2,:), color = 'r',
        sol2(1,:),sol2(2,:), color='g', ptoTg(1),ptoTg(2),color='g*');
  xlim([0 3])
  ylim([0 5])
  grid on
  pause(0.01)
endfor
