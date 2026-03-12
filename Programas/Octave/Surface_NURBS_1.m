% --- CARGA DE DATOS DE LA SUPERFICIE ---
dataX = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj3RefinedX.txt', ' ');
dataY = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj3RefinedY.txt', ' ');
dataZ = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Superficies\SuperficieRationalEj3RefinedZ.txt', ' ');

% --- CARGA DE DATOS DE LOS PUNTOS DE CONTROL ---
dataCtrlX = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej3RefinedX.txt', ' ');
dataCtrlY = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej3RefinedY.txt', ' ');
dataCtrlZ = dlmread('C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\CtrlPts\Ej3RefinedZ.txt', ' ');

[n1, n2] = size(dataCtrlX);

% --- CONFIGURACIÓN DE LA FIGURA ---
figure('Position', [100, 100, 800, 600]);
hold on; grid on;

% Trazado de la superficie NURBS
surf(dataX, dataY, dataZ);
% colormap(parula); % Descomenta esto si prefieres un mapa de colores más moderno que el jet por defecto

% Trazado de la red de control (líneas rojas discontinuas con marcadores)
for i = 1:n1
  plot3(dataCtrlX(i,:), dataCtrlY(i,:), dataCtrlZ(i,:), 'r--o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'r');
end
for j = 1:n2
  plot3(dataCtrlX(:,j), dataCtrlY(:,j), dataCtrlZ(:,j), 'r--o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'r');
end

% --- AJUSTE DINÁMICO DE EJES ---
min_X = floor(min([dataX(:); dataCtrlX(:)]));
max_X = ceil(max([dataX(:); dataCtrlX(:)]));
min_Y = floor(min([dataY(:); dataCtrlY(:)]));
max_Y = ceil(max([dataY(:); dataCtrlY(:)]));
min_Z = floor(min([dataZ(:); dataCtrlZ(:)]));
max_Z = ceil(max([dataZ(:); dataCtrlZ(:)]));

xlim([min_X, max_X]);
ylim([min_Y, max_Y]);
zlim([min_Z, max_Z]);

% Generar 5 marcas equidistantes para cada eje
set(gca, 'XTick', linspace(min_X, max_X, 5));
y_ticks = linspace(min_Y, max_Y, 5);
set(gca, 'YTick', y_ticks);
set(gca, 'ZTick', linspace(min_Z, max_Z, 5));

set(gca, 'TickDir', 'out', 'TickLength', [0.1, 0.1]);

% Truco TeX para desplazar los números del eje Y
y_labels = cell(length(y_ticks), 1);
for i = 1:length(y_ticks)
  y_labels{i} = sprintf('%g~~', y_ticks(i));
end
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'YTickLabel', y_labels);

xlabel(sprintf('X\n\n\n'));
ylabel(sprintf('\n\n\n Y'));
zlabel(sprintf('Z\n\n\n'));

view(3);
set(gca, 'Position', [0.15 0.15 0.7 0.8]);

hold off;

% --- EXPORTACIÓN A SVG ---
print('Superficie_NURBS_Ej2_Refined.svg', '-dsvg');
disp('Superficie exportada exitosamente como Superficie_NURBS_Ej1.svg');
