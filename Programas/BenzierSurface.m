%BenzierSurface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constantes de control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=100; %numero de pasos
u=0:1/m:1; %vector del parametro de la curva

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINIMOS LOS P (puntos de control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=zeros(3, 4, 5);

P(:,:,1)=[0, 0, 0, 0;
          0, 1, 2, 3;
          0, 2, 3, 3];

P(:,:,2)=[1, 1, 1, 1;
          0, 1, 2, 3;
          0, 2, 3, 3];

P(:,:,3)=[2, 2, 2, 2;
          0, 3, 4, 4;
          1, 2, 3, 3];

P(:,:,4)=[3, 3, 3, 3;
          0, 3, 3, 5;
          1, 2, 3, 3];

P(:,:,5)=[4, 4, 4, 4;
          0, 3, 4, 5;
          5, 5, 5, 5];


[aux n1 n2] = size(P); %guardamos el grado n1 (tener en cuenta que sera n1+1 realmente
                       %guardamos el numero de particiones que hacemos sobre el eje x n2
                       % ya que octave indexiza desde el 1)
sol=zeros(aux,m+1,m+1); %sera nuestra superficie
sol2=zeros(aux,m+1);%sera el vector plano tangente

for k1=1:m+1    %iterador de Vo
  for k2=1:m+1  %iterador de Uo
    Bold1=zeros(n1, 1);  %reinicio de las variables para cada punto de la curva
    Bold1(1)=1; %Caso B(0,0) n=0
    Bnew1=zeros(n1,1);

    for j1=2:n1
      Bnew1=zeros(n1,1);
      for i=1:j1-1
        Bnew1(i)=Bnew1(i) + (1-u(k1))*Bold1(i);
        Bnew1(i+1)=Bnew1(i+1) + (u(k1))*Bold1(i);
      endfor
      Bold1=Bnew1;
    endfor

    Bold2=zeros(n2, 1);
    Bold2(1)=1;
    Bnew2=zeros(n2,1);


    for j2=2:n2
      Bnew2=zeros(n2,1);
      for i=1:j2-1
        Bnew2(i)=Bnew2(i) + (1-u(k2))*Bold2(i);
        Bnew2(i+1)=Bnew2(i+1) + (u(k2))*Bold2(i);
      endfor
      Bold2=Bnew2;
    endfor


    for j1=1:n1
      for j2=1:n2
        sol(:,k1, k2)=sol(:,k1, k2)+Bnew1(j1)*Bnew2(j2)*P(:,j1,j2);
      endfor
    endfor


  endfor
endfor


  for i=1:m+1
    %for uo=0:0.01:1 %bucle que mueve la tangente y el punto tg

    %    rectaTg=P;
    %    for j=1:n-2
    %      for i=1:n-j
    %        rectaTg(:,i)=(1-uo)*rectaTg(:,i)+(uo)*rectaTg(:,i+1);
    %      endfor
    %    endfor
    %    ptoTg=rectaTg(:,1)*(1-uo)+rectaTg(:,2)*(uo);

      plot3(sol(1,:,i),sol(2,:,i), sol(3,:,i), color='b');
      hold on
            %rectaTg(1,1:2),rectaTg(2,1:2), rectaTg(3,1:2), color='g',
            %ptoTg(1), ptoTg(2), ptoTg(3), color='g*');
      xlim([0 5])
      ylim([0 5])
      zlim([0 5])
      grid on

    %endfor

  endfor

  for j=1:n2
    plot3(P(1,:,j),P(2,:,j), P(3,:,j), color='r');
  endfor
  %fix para el plot de octave
  Paux=zeros(aux,n2,n1);
  for(k=1:aux)
    for i=1:n1
      for j=1:n2
        Paux(k,j,i)=P(k,i,j);
      endfor
    endfor
  endfor

  for i=1:n1
   plot3(Paux(1,:,i),Paux(2,:,i), Paux(3,:,i), color='r');
  endfor

  #rectaTg1=P;
  #for j=1:n2-2
  #  for i=1:n2-j
  #    rectaTg1(:,:,i)=(1-uo)*rectaTg1(:,:,i)+(uo)*rectaTg1(:,:,i+1);
  #  endfor
  #endfor
  #ptoTg1=rectaTg1(:,:,1)*(1-uo)+rectaTg1(:,:,2)*(uo);
  #
  #P3=plot3(rectaTg1(1,:,i),rectaTg1(2,:,i), rectaTg1(3,:,i), color='g');
  #G+P1+P2
  #pause(0.05)
  #hold off

#plot3(rectaTg1(1,:,i),rectaTg1(2,:,i), rectaTg1(3,:,i), color='g');
hold off
