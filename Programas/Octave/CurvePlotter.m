% ##########################################################################
% SCRIPT OCTAVE PARA TRAZAR CURVA NURBS 3D Y PUNTOS DE CONTROL
% Versión FINAL: Asume que la data está en 3 filas x N columnas (3 x 1001)
% ##########################################################################

% --- CONFIGURACIÓN DE RUTAS Y ARCHIVOS ---
RUTA_ARCHIVO_CURVA = 'C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/EjCurva3d.txt';
NOMBRE_ARCHIVO_SVG = 'CurvaNurbs_3D_Plot.svg';

% 1. Definir la matriz de Puntos de Control (PC)
PC = [
    10, 0, 0;    % P0
    15, 10, 5;   % P1
    5, 15, 20;   % P2 (Punto con peso alto)
    0, 5, 10;    % P3
    10, 0, 0     % P4
];

disp('Cargando datos de la curva...');

% 2. Leer los datos de la curva desde el archivo
try
    % D será una matriz de M x K. En su caso, D será 3 x 1001.
    D = load(RUTA_ARCHIVO_CURVA);
catch
    error('ERROR: No se pudo cargar el archivo. Comprueba que la ruta sea correcta: %s', RUTA_ARCHIVO_CURVA);
end

% 3. Procesar los datos: Transposición y verificación
[M, K] = size(D);

if M == 3 && K >= 2
    % CASO DEFINITIVO: 3 filas x N columnas (e.g., 3 x 1001).
    % Las coordenadas (X, Y, Z) están en las filas, y los puntos en las columnas.
    C = D'; % Transponer la matriz D. C ahora es N x 3 (1001 x 3).
    N = K; % El número de puntos es el número de columnas originales (1001).
    disp(['Formato detectado: [X; Y; Z] en 3 filas x ', num2str(N), ' columnas. Transpuesto a N x 3.']);

elseif K == 3
    % Caso Secundario: Data en N filas x 3 columnas (N x 3).
    C = D;
    N = M;
    disp(['Formato detectado: [X Y Z] en ', num2str(N), ' filas x 3 columnas.']);

else
    % Caso de error: Número de filas/columnas no compatible (no es 3xN ni Nx3)
    error('ERROR DE FORMATO: Se esperaba una matriz 3xN o Nx3. Se detectó una matriz %d x %d.', M, K);
end

% 4. Configuración de la figura y el trazado
figure;
hold on;
grid on;

% Trazar la curva NURBS (azul)
if N > 1
    plot3(C(:,1), C(:,2), C(:,3), 'b-', 'LineWidth', 2);
    disp('Curva trazada exitosamente.');
else
    disp('Advertencia: El archivo de curva no contiene suficientes puntos (N < 2) para dibujar una línea.');
end

% Trazar el Polígono de Control (líneas rojas discontinuas)
plot3(PC(:,1), PC(:,2), PC(:,3), 'r--', 'LineWidth', 1, 'DisplayName', 'Polígono de Control');

% Trazar los Puntos de Control (círculos rojos)
plot3(PC(:,1), PC(:,2), PC(:,3), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', 'Puntos de Control');

% 5. Configuración de la visualización
%title('Curva NURBS 3D con Polígono de Control', 'FontSize', 14);
%xlabel('Eje X', 'FontSize', 12);
%ylabel('Eje Y', 'FontSize', 12);
%zlabel('Eje Z', 'FontSize', 12);

% Ajustar límites de los ejes (usando la matriz C transpuesta correctamente)
min_coords = min([C; PC]);
max_coords = max([C; PC]);
margen = 2;
axis([min_coords(1)-margen max_coords(1)+margen ...
      min_coords(2)-margen max_coords(2)+margen ...
      min_coords(3)-margen max_coords(3)+margen]);

view(3);

set(gca, 'Box', 'on');

% 6. Guardar en formato SVG
print(NOMBRE_ARCHIVO_SVG, '-dsvg');
disp(['Gráfico final guardado en: ', fullfile(pwd, NOMBRE_ARCHIVO_SVG)]);
