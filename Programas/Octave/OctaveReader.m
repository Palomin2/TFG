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
  data1=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\RatBasisFuncts.txt');
  data2=dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\RatBasisFuncts_Derivates.txt');

  [n1 n2] = size(data2);
  hold on
  for i=1:n1
    plot(data1(i,:), color='r');
    #pause(0.5)
  endfor
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
  data1 = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\BasisFuncts2.txt');

  [n1 n2] = size(data1);
  colores = {'r', 'g', 'b', 'm', 'c', [1, 0.5, 0], [0.5, 0.5, 0.5], [0.2, 0.8, 0.2]};
  num_colores = length(colores);

  figure;
  % --- AUMENTAMOS LOS MÁRGENES DE LA FIGURA ---
  % Dejamos el mismo margen del 15% por la izquierda y por abajo
  set(gca, 'Position', [0.15 0.15 0.8 0.8]);
  % --------------------------------------------

  hold on;
  grid on;

  % --- AÑADIMOS ESPACIO A LOS LADOS (10% del ancho) ---
  margen_x = 0.1 * (n2 - 1);
  limite_izq = 1 - margen_x;
  limite_der = n2 + margen_x;

  axis([limite_izq limite_der -0.2 1.2]);
  box on;
  % ----------------------------------------------------

  % --- BLOQUE PARA DESPLAZAR LOS NÚMEROS DE LOS EJES ---
  posiciones_ticks_fisicas = [1, 21, 41, 61, 81, 101];
  etiquetas_deseadas = {'0', '0.2', '0.4', '0.6', '0.8', '1'};
  yticks_vals = [0, 0.2, 0.4, 0.6, 0.8, 1];

  set(gca, 'XTick', posiciones_ticks_fisicas);
  set(gca, 'YTick', yticks_vals);
  set(gca, 'XTickLabel', []);
  set(gca, 'YTickLabel', []);

  % Ajustamos el desplazamiento:
  desplazamiento_y = -0.32;
  desplazamiento_x = limite_izq - 10; % Lo alejamos 9 puntos desde el nuevo borde izquierdo

  for k = 1:length(posiciones_ticks_fisicas)
      text(posiciones_ticks_fisicas(k), desplazamiento_y, etiquetas_deseadas{k}, ...
          'HorizontalAlignment', 'center', 'Clipping', 'off');
  end

  for k = 1:length(yticks_vals)
      text(desplazamiento_x, yticks_vals(k), sprintf('%g', yticks_vals(k)), ...
          'HorizontalAlignment', 'right', 'Clipping', 'off');
  end
  % -----------------------------------------------------

  for i = 1:n1
      color_index = mod(i - 1, num_colores) + 1;
      plot(data1(i,:), 'color', colores{color_index}, 'LineWidth', 1.5);
  end
  hold off;



  nombre_archivo_svg = 'funciones_base_NURBS_2.svg';
  print(nombre_archivo_svg, '-dsvg');
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

  scatter3(dataCtrlX(:), dataCtrlY(:), dataCtrlZ(:), 30, 'k', 'filled', 'MarkerEdgeColor', 'w');
  view(30, 45);
  grid on;

  print -dsvg 'Superficie_NURBS_Ej1.svg';
  hold off;
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
  scatter3(dataCtrlX(:), dataCtrlY(:), dataCtrlZ(:), 30, 'k', 'filled', 'MarkerEdgeColor', 'w');
  view(30, 45);
  grid on;

  print -dsvg 'Superficie_NURBS_Ej2.svg';
  hold off;
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
  ;

  for i=1:n1
    plot3(dataCtrlX(i,:),dataCtrlY(i,:),dataCtrlZ(i,:), color='r');
  endfor

  for j=1:n2
    plot3(dataCtrlX(:,j),dataCtrlY(:,j),dataCtrlZ(:,j), color='r');
  endfor
  scatter3(dataCtrlX(:), dataCtrlY(:), dataCtrlZ(:), 30, 'k', 'filled', 'MarkerEdgeColor', 'w');
  view(30, 45);
  grid on;

  print -dsvg 'Superficie_NURBS_Ej2_Refined.svg';
  hold off
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
  data1 = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConst_h=16_p=4_test1.txt');
  data2 = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConst_Analytic_test1.txt');
  data3 = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConstDers_h=16_p=4_test1.txt');
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
  print(h1,'ej1_p4_h16_sol.svg','-dsvg');

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
  print(h2,'ej1_p4_h16_ders.svg','-dsvg');
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

   elseif(it==19)
   figure
  dataX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEjSuperficieTrivX.txt', ' ');
  dataY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEjSuperficieTrivY.txt', ' ');
  dataZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEjSuperficieTrivZ.txt', ' ');

  hold on
  surf(dataX, dataY, dataZ);

  dataCtrlX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\EjSuperficieTrivX.txt', ' ');
  dataCtrlY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\EjSuperficieTrivY.txt', ' ');
  dataCtrlZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\EjSuperficieTrivZ.txt', ' ');

  [n1 n2] = size(dataCtrlX);

  for i=1:n1
    plot3(dataCtrlX(i,:),dataCtrlY(i,:),dataCtrlZ(i,:), color='r');
  endfor

  for j=1:n2
    plot3(dataCtrlX(:,j),dataCtrlY(:,j),dataCtrlZ(:,j), color='r');
  endfor
  figure
  dataX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEjSuperficieTrivRefinementX.txt', ' ');
  dataY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEjSuperficieTrivRefinementY.txt', ' ');
  dataZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEjSuperficieTrivRefinementZ.txt', ' ');

  hold on
  surf(dataX, dataY, dataZ);

  dataCtrlX =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\EjSuperficieTrivRefinementX.txt', ' ');
  dataCtrlY =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\EjSuperficieTrivRefinementY.txt', ' ');
  dataCtrlZ =  dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\EjSuperficieTrivRefinementZ.txt', ' ');

  [n1 n2] = size(dataCtrlX);

  for i=1:n1
    plot3(dataCtrlX(i,:),dataCtrlY(i,:),dataCtrlZ(i,:), color='r');
  endfor

  for j=1:n2
    plot3(dataCtrlX(:,j),dataCtrlY(:,j),dataCtrlZ(:,j), color='r');
  endfor
  elseif (it == 20)
  %% --- Parámetros configurables ---
  h_val = 32;    % resolución
  p_val = 2;       % grado p

  base_path = 'C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\';
  test_name = 'testCircle';

  %% --- Cargar y preparar datos ---
  data1 = dlmread([base_path sprintf('EjSinNonConst_h=%d_p=%d_%s.txt', h_val, p_val, test_name)]);
  data2 = dlmread([base_path sprintf('EjSinNonConst_Analytic_%s.txt', test_name)]);
  data3 = dlmread([base_path sprintf('EjSinNonConstDers_h=%d_p=%d_%s.txt', h_val, p_val, test_name)]);
  data4 = dlmread([base_path sprintf('EjSinNonConstDers_Analytic_%s.txt', test_name)]);

  Length = 1;
  data1 = data1 / (Length * Length);
  data3 = data3 * pi;

  [n1, n2] = size(data1);

  [a, b] = deal(0, 1);  % intervalo
  [xq, wq] = gauss_legendre(n2, a, b);
  %wq = wq / 2;

  %% --- Figura 1: Solución ---
  h1 = figure('Position', [100 100 800 600]);
  hold on
  for i = 1:n1
      plot(xq, data1(i,:), 'r');
  endfor
  for i = 1:n1
      plot(xq, data2(i,:), 'b');
  endfor
  legend({'$u_{\text{num}}$', '$u_{\text{ana}}$'}, 'Interpreter', 'latex', 'FontSize', 30);
  title('Solución', 'FontSize', 30, 'Interpreter', 'tex');
  xlabel('x', 'FontSize', 30, 'Interpreter', 'latex');
  ylabel('u(x)', 'FontSize', 30, 'Interpreter', 'latex');

  set(gca, 'LooseInset', get(gca, 'TightInset'));
  print(h1, sprintf('ejCircle_p%d_h%d_sol.svg', p_val, h_val), '-dsvg');
  hold off

  %% --- Figura 2: Derivada ---
  h2 = figure('Position', [100 100 800 600]);
  hold on
  for i = 1:n1
      plot(xq, data3(i,:), 'g');
  endfor
  for i = 1:n1
      plot(xq, data4(i,:), 'y');
  endfor
  legend({'$u^{\prime}_{\text{num}}$', '$u^{\prime}_{\text{ana}}$'}, 'Interpreter', 'latex', 'FontSize', 30);
  title('Derivada', 'FontSize', 30, 'Interpreter', 'tex');
  xlabel('x', 'FontSize', 30, 'Interpreter', 'latex');
  ylabel('$u^{\prime}(x)$', 'FontSize', 30, 'Interpreter', 'latex');

  set(gca, 'LooseInset', get(gca, 'TightInset'));
  print(h2, sprintf('ejCircle_p%d_h%d_ders.svg', p_val, h_val), '-dsvg');
  hold off

  %% --- Cálculo de normas ---
  [a, b] = deal(0, 1);
  [xq, wq] = gauss_legendre(n2, a, b);
  %xq = (xq + 1);
  L2sq = sum(sum((data1 - data2).^2 .* wq));
  L2norm = sqrt(L2sq)
  H1sqs = sum(sum((data3 - data4).^2 .* wq));
  H1norm = sqrt(L2sq + H1sqs)
  normInf = norm(data1 - data2, Inf)
  %% --- Test Cuadratura ---
  %[a,b]=deal(0,2);
  %[xq,wq]=gauss_legendre(n2,a,b);
  %I1 = sum(wq);     % debe ser 2
  %Ix = sum(wq .* xq); % debe ser integral x dx over [0,2] = 2
  %fprintf('sum(wq)=%.12g (esperado 2), sum(wq.*xq)=%.12g (esperado 2)\n', I1, Ix);

endif
hold off
