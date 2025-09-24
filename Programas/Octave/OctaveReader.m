function OctaveReader(it)
#it=2
%figure
%h=figure('Position',[100 100 600 400]);  % tamaño de ventana: ancho x alto en pixeles
if(it==1)
  dataX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\BezierSurf\TestFicheroDatosBezierSurfaceX.txt', ' ');
  dataY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\BezierSurf\TestFicheroDatosBezierSurfaceY.txt', ' ');
  dataZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\BezierSurf\TestFicheroDatosBezierSurfaceZ.txt', ' ');

  hold on
  surf(dataX, dataY, dataZ);

  dataCtrlX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\BezierSurf\TestFicheroDatosPtosCtrlX.txt', ' ');
  dataCtrlY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\BezierSurf\TestFicheroDatosPtosCtrlY.txt', ' ');
  dataCtrlZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\BezierSurf\TestFicheroDatosPtosCtrlZ.txt', ' ');

  [n1 n2] = size(dataCtrlX);

  for i=1:n1
    plot3(dataCtrlX(i,:),dataCtrlY(i,:),dataCtrlZ(i,:), color='r');
  endfor

  for j=1:n2
    plot3(dataCtrlX(:,j),dataCtrlY(:,j),dataCtrlZ(:,j), color='r');
  endfor
elseif(it==2)
  data1=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\BasisFuncts.txt');
  data2=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\RatBasisFuncts_Derivates.txt');

  [n1 n2] = size(data2);
  hold on
  for i=1:n1
    plot(data2(i,:), color='r');
    #pause(0.5)
  endfor

  plot(data2(4,:), color='b');
  hold off
elseif(it==3)
  data1=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Curvas\CurvaEj1.txt');
  [n1 n2] = size(data1);
  data2=[0, 1, 1, 2, 3, 2.5, 1, 1, 2, 3, 2, 3, 2, 1, 2;
         0, 0, 1, 1, 0, -1, -0.5, 0, 1, 2, -1, -1, 0, 0, 1];
  [n3 n4] = size(data2);
  hold on

  for i=1:n2
    plot(data1(1,:),data1(2,:));
    #pause(0.5)
  endfor
  for i=1:n4
    plot(data2(1,:),data2(2,:), color='b');
    plot(data2(1,:),data2(2,:), color='b*');
    #pause(0.5)
  endfor
  hold off
elseif(it==9)
  data3=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Curvas\CurvaEj1DersAlg1.txt');
  data4=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Curvas\CurvaEj1DersAlg2.txt');
  [n1 n2] = size(data3);
    for i=1:n2
    plot(data3(1,:),data3(2,:),color='y');
    #pause(0.5)
  endfor

  [n1 n2] = size(data4);
  for i=1:n2
    #plot(data4(1,:),data4(2,:),color='g');
    #pause(0.5)
  endfor



  hold off
elseif(it==4)
  dataX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieEj2X.txt', ' ');
  dataY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieEj2Y.txt', ' ');
  dataZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieEj2Z.txt', ' ');

  hold on
  surf(dataX, dataY, dataZ);

  dataCtrlX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej2X.txt', ' ');
  dataCtrlY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej2Y.txt', ' ');
  dataCtrlZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej2Z.txt', ' ');

  [n1 n2] = size(dataCtrlX);

  for i=1:n1
    plot3(dataCtrlX(i,:),dataCtrlY(i,:),dataCtrlZ(i,:), color='r');
  endfor

  for j=1:n2
    plot3(dataCtrlX(:,j),dataCtrlY(:,j),dataCtrlZ(:,j), color='r');
  endfor
elseif(it==5)
  data1=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\AllBasisFuncts.txt');
  %data2=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\BasisFuncts_Derivates.txt');

  [n1 n2] = size(data1);
  hold on
  for i=1:n1
    plot(data1(i,:), color='r');
    #pause(0.5)
  endfor

  %plot(data1(4,:), color='b');
  %hold off
elseif(it==6)
  data1=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Curvas\CurvaEj1Nurbs.txt');
  [n1 n2] = size(data1);
  data2=[0, 1, 1, 2, 3, 2.5, 1;
         0, 0, 1, 1, 0, -1, -0.5];
  [n3 n4] = size(data2);
  hold on
  for i=1:n2
  plot(data1(1,:),data1(2,:));
  #pause(0.5)
  endfor
  for i=1:n4
    plot(data2(1,:),data2(2,:), color='b');
    plot(data2(1,:),data2(2,:), color='b*');
    #pause(0.5)
  endfor

  elseif(it==7)
  data1=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\RationalBasisFuncts.txt');


  [n1 n2] = size(data1);
  hold on
  for i=1:n1
    plot(data1(i,:), color='r');
    #pause(0.5)
  plot(data1(2,:), color='b');
  endfor
  hold off
  elseif(it==8)
  dataX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj2X.txt', ' ');
  dataY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj2Y.txt', ' ');
  dataZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj2Z.txt', ' ');

  hold on
  surf(dataX, dataY, dataZ);

  dataCtrlX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej2X.txt', ' ');
  dataCtrlY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej2Y.txt', ' ');
  dataCtrlZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej2Z.txt', ' ');

  [n1 n2] = size(dataCtrlX);

  for i=1:n1
    plot3(dataCtrlX(i,:),dataCtrlY(i,:),dataCtrlZ(i,:), color='r');
  endfor

  for j=1:n2
    plot3(dataCtrlX(:,j),dataCtrlY(:,j),dataCtrlZ(:,j), color='r');
  endfor
  elseif(it==10)
  data1=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Curvas\CurvaEj1CornerCut.txt');
  [n1 n2] = size(data1);
  data2=[0, 1, 3, 2, 4, 3, 5;
         0, 1, 1, -1, -1, 0, 0];
  [n3 n4] = size(data2);
  hold on
  for i=1:n2
  plot(data1(1,:),data1(2,:));
  #pause(0.5)
  endfor
  for i=1:n4
    #plot(data2(1,:),data2(2,:), color='b');
    #plot(data2(1,:),data2(2,:), color='b*');
    #pause(0.5)
  endfor
  elseif(it==11)
  dataX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj2InsertionX.txt', ' ');
  dataY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj2InsertionY.txt', ' ');
  dataZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj2InsertionZ.txt', ' ');

  hold on
  surf(dataX, dataY, dataZ);

  dataCtrlX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej2InsertionX.txt', ' ');
  dataCtrlY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej2InsertionY.txt', ' ');
  dataCtrlZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej2InsertionZ.txt', ' ');

  [n1 n2] = size(dataCtrlX);

  for i=1:n1
    plot3(dataCtrlX(i,:),dataCtrlY(i,:),dataCtrlZ(i,:), color='r');
  endfor

  for j=1:n2
    plot3(dataCtrlX(:,j),dataCtrlY(:,j),dataCtrlZ(:,j), color='r');
  endfor
  elseif(it==12)
  data1=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Curvas\CurvaEj1RefineKnotVect.txt');
  [n1 n2] = size(data1);
  data2=[0.000,  0.500,  0.875,  0.625,  0.875,  1.625,  2.250,  2.750,  2.875,  2.625,  1.750,  1.000;
0.000  0.000,  0.125,  0.375,  0.625,  0.875,  0.750,  0.250,  -0.250,  -0.750,  -0.750,  -0.500];
  [n3 n4] = size(data2);
  hold on
  for i=1:n2
  plot(data1(1,:),data1(2,:));
  #pause(0.5)
  endfor
  for i=1:n4
    plot(data2(1,:),data2(2,:), color='b');
    plot(data2(1,:),data2(2,:), color='b*');
    #pause(0.5)
  endfor
    elseif(it==13)
  dataX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj2RefineKnotVectSurfX.txt', ' ');
  dataY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj2RefineKnotVectSurfY.txt', ' ');
  dataZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj2RefineKnotVectSurfZ.txt', ' ');

  hold on
  surf(dataX, dataY, dataZ);

  dataCtrlX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\SuperficieRationalEj2RefineKnotVectSurfX.txt', ' ');
  dataCtrlY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\SuperficieRationalEj2RefineKnotVectSurfY.txt', ' ');
  dataCtrlZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\SuperficieRationalEj2RefineKnotVectSurfZ.txt', ' ');

  [n1 n2] = size(dataCtrlX);

  for i=1:n1
    plot3(dataCtrlX(i,:),dataCtrlY(i,:),dataCtrlZ(i,:), color='r');
  endfor

  for j=1:n2
    plot3(dataCtrlX(:,j),dataCtrlY(:,j),dataCtrlZ(:,j), color='r');
  endfor
  elseif(it==14)
  dataX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj3X.txt', ' ');
  dataY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj3Y.txt', ' ');
  dataZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj3Z.txt', ' ');

  hold on
  surf(dataX, dataY, dataZ);

  dataCtrlX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej3X.txt', ' ');
  dataCtrlY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej3Y.txt', ' ');
  dataCtrlZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej3Z.txt', ' ');

  [n1 n2] = size(dataCtrlX);

  for i=1:n1
    plot3(dataCtrlX(i,:),dataCtrlY(i,:),dataCtrlZ(i,:), color='r');
  endfor

  for j=1:n2
    plot3(dataCtrlX(:,j),dataCtrlY(:,j),dataCtrlZ(:,j), color='r');
  endfor
  elseif(it==15)
  dataX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj3RefinedX.txt', ' ');
  dataY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj3RefinedY.txt', ' ');
  dataZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj3RefinedZ.txt', ' ');

  hold on
  surf(dataX, dataY, dataZ);

  dataCtrlX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej3RefinedX.txt', ' ');
  dataCtrlY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej3RefinedY.txt', ' ');
  dataCtrlZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej3RefinedZ.txt', ' ');

  [n1 n2] = size(dataCtrlX);

  for i=1:n1
    plot3(dataCtrlX(i,:),dataCtrlY(i,:),dataCtrlZ(i,:), color='r');
  endfor

  for j=1:n2
    plot3(dataCtrlX(:,j),dataCtrlY(:,j),dataCtrlZ(:,j), color='r');
  endfor
  elseif(it==16)
  data1=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSin_h=80_p=2.txt');
  data2=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSin_Analytic_h=80.txt');
  data3=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinDers_h=80_p=2.txt');
  data4=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinDers_Analytic_h=80.txt');
  [n1 n2] = size(data1);
  hold on
  for i=1:n1
    plot(data1(i,:), color='r');
    #pause(0.5)
  endfor
  for i=1:n1
    plot(data2(i,:), color='b');
    #pause(0.5)
  endfor
  hold off               % número de puntos de Gauss
  [a,b] = deal(0,1);   % intervalo
  [xq, wq] = gauss_legendre(n2, a, b);
  wq=wq/2;
  L2sq = sum( wq .* (data1 - data2).^2 );
  L2norm= sqrt(L2sq)
  H1sqs = sum( wq .* (data3 - data4).^2 );
  H1norm = sqrt(L2sq + H1sqs)
  normInf = norm(data1 - data2, Inf)
  elseif(it==17)
  data1=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConst_h=40_p=4.txt');
  data2=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConst_Analytic_h=40.txt');
  data3=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConstDers_h=40_p=4.txt');
  data4=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConstDers_Analytic_h=40.txt');
  [n1 n2] = size(data1);
  hold on
  for i=1:n1
    plot(data1(i,:), color='r');
    #pause(0.5)
  endfor
  for i=1:n1
    plot(data2(i,:), color='b');
    #pause(0.5)
  endfor
  figure
  hold on
    for i=1:n1
    plot(data3(i,:), color='g');
    #pause(0.5)
  endfor
  for i=1:n1
    plot(data4(i,:), color='y');
    #pause(0.5)
  endfor
  hold off               % número de puntos de Gauss
  [a,b] = deal(0,1);   % intervalo
  [xq, wq] = gauss_legendre(n2, a, b);
  xq = (xq + 1) / 2
  wq=wq/2;
  L2sq = sum( wq .* (data1 - data2).^2 );
  L2norm= sqrt(L2sq)
  H1sqs = sum( wq .* (data3 - data4).^2 );
  H1norm = sqrt(L2sq + H1sqs)
  normInf = norm(data1 - data2, Inf)
  elseif(it==18)
  %% --- Cargar y preparar datos ---
  data1 = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConst_h=32_p=4_test1.txt');
  data2 = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConst_Analytic_test1.txt');
  data3 = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConstDers_h=32_p=4_test1.txt');
  data4 = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConstDers_Analytic_test1.txt');

  Length = 1;
  data1 = data1 / (Length*Length);
  data3 = data3 / Length;

  [n1, n2] = size(data1);

  [a,b] = deal(0,1);  % intervalo
  [xq, wq] = gauss_legendre(n2, a, b);
  wq = wq / 2;

  %% --- Figura 1: Solución ---
  h1 = figure('Position',[100 100 800 600]);
  hold on
  for i=1:n1
      plot(xq, data1(i,:), 'r');
  endfor
  for i=1:n1
      plot(xq, data2(i,:), 'b');
  endfor
  legend({'$u_{\text{num}}$','$u_{\text{ana}}$'}, 'Interpreter','latex','FontSize',30);
  title('Solución', 'FontSize', 30, 'Interpreter', 'tex');
  xlabel('x', 'FontSize', 30, 'Interpreter', 'latex');
  ylabel('u(x)', 'FontSize', 30, 'Interpreter', 'latex');

  % --- Ajustar márgenes tight para SVG ---
  set(gca,'LooseInset',get(gca,'TightInset'));
  % Exportar en SVG
  print(h1,'ej1_p4_h32_sol.svg','-dsvg');

  hold off

  %% --- Figura 2: Derivada ---
  h2 = figure('Position',[100 100 800 600]);
  hold on
  for i=1:n1
      plot(xq, data3(i,:), 'g');
  endfor
  for i=1:n1
      plot(xq, data4(i,:), 'y');
  endfor
  legend({'$u^{\prime}_{\text{num}}$','$u^{\prime}_{\text{ana}}$'}, 'Interpreter','latex','FontSize',30);
  title('Derivada', 'FontSize', 30, 'Interpreter', 'tex');
  xlabel('x', 'FontSize', 30, 'Interpreter', 'latex');
  ylabel('$u^{\prime}(x)$', 'FontSize', 30, 'Interpreter', 'latex');

  set(gca,'LooseInset',get(gca,'TightInset'));
  print(h2,'ej1_p4_h32_ders.svg','-dsvg');
  hold off

  %% --- Figura 3: Curva coloreada ---
  %curva   = load('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Curvas\CurvaEjPruebaFEM2.txt');
  %valores = data1;

  %x = curva(1,:);
  %y = curva(2,:);
  %z = curva(3,:);
  %h3 = figure('Position',[100 100 800 600]);
  %scatter3(x, y,z, 60, valores, 'filled');
  %colormap(jet);
  %colorbar;
  %axis equal;
  %grid on;

  %xlabel('x', 'FontSize',14);
  %ylabel('y', 'FontSize',14);
  %zlabel('z', 'FontSize',14);

 % Ajuste del título para que no se corte (Octave)
  %t = title('Curva en $R^3$', 'Interpreter','latex','FontSize',30);
  %set(gca,'LooseInset',[0.1 0.1 0.1 0.25]);  % margen superior más grande

  %pos = get(t, 'position');   % vector [x y z]
  %pos(2) = pos(2) - 0.05;     % subir el título ligeramente
  %set(t, 'position', pos);


  % Exportar a SVG
  %print(h3,'ej1_p4_h8_curve.svg','-dsvg');
  data1
  [a,b] = deal(0,2);   % intervalo
  [xq, wq] = gauss_legendre(n2, a, b);
  xq = (xq + 1);
  wq=wq;
  L2sq = sum( wq .* (data1 - data2).^2 );
  L2norm= sqrt(L2sq)
  H1sqs = sum( wq .* (data3 - data4).^2 );
  H1norm = sqrt(L2sq + H1sqs)
  normInf = norm(data1 - data2, Inf)
endif
hold off
