% ##########################################################################
% SCRIPT OCTAVE PARA TRAZAR CURVA NURBS 3D Y PUNTOS DE CONTROL
% Versión FINAL: Asume que la data está en 3 filas x N columnas (3 x 1001)
% ##########################################################################

% --- CONFIGURACIÓN DE RUTAS Y ARCHIVOS ---
RUTA_ARCHIVO_CURVA = 'C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/EjCurvaCirculo2.txt';
NOMBRE_ARCHIVO_SVG = 'CurvaCirculoNURBS_3D_Plot.svg';

% 1. Definir la matriz de Puntos de Control (PC)
PC = [
    1, 0, 0;    % P0
    1, 1, 0;   % P1
    0, 1, 0;   % P2 (Punto con peso alto)
    -1, 1, 0;    % P3
    -1, 0, 0;     % P4
    -1, -1, 0;
    0, -1, 0;
    1, -1, 0;
    1, 0, 0
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
figure('Position', [100, 100, 800, 600]);
hold on;
grid on;

% Trazar la curva NURBS (azul)
if N > 1
    plot3(C(:,1), C(:,2), C(:,3), 'b-', 'LineWidth', 2);
    disp('Curva trazada exitosamente.');
else
    disp('Advertencia: El archivo de curva no contiene suficientes puntos.');
end

% Trazar el Polígono y Puntos de Control (rojo discontinuo con círculos)
plot3(PC(:,1), PC(:,2), PC(:,3), 'r--o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'r');

% 5. Configuración de la visualización (EJES DINÁMICOS Y LIMPIOS)
min_coords = min([C; PC]); % Vector 1x3 con los mínimos [X, Y, Z]
max_coords = max([C; PC]); % Vector 1x3 con los máximos [X, Y, Z]

% Redondear a enteros para tener límites de caja perfectos
min_X = floor(min_coords(1));
max_X = ceil(max_coords(1));
min_Y = floor(min_coords(2));
max_Y = ceil(max_coords(2));
min_Z = floor(min_coords(3));
max_Z = ceil(max_coords(3));

% CONDICIONAL DE SEGURIDAD PARA CURVAS PLANAS
if min_Z == max_Z
    min_Z = min_Z - 1;
    max_Z = max_Z + 1;
end

% Establecer límites cerrados a números enteros
xlim([min_X, max_X]);
ylim([min_Y, max_Y]);
zlim([min_Z, max_Z]);

% Forzar las marcas de los ejes a saltos de 5 en 5 (dado que el rango es ~20)
set(gca, 'XTick', min_X:1:max_X);
y_ticks = min_Y:1:max_Y;
set(gca, 'YTick', y_ticks);
z_ticks = linspace(min_Z, max_Z, 3);
set(gca, 'ZTick', z_ticks);

% Sacar las marcas hacia afuera con una longitud prudente
set(gca, 'TickDir', 'out', 'TickLength', [0.1, 0.1]);

% ** TRUCO CON TeX: Espacios insecables a la derecha para el eje Y **
y_labels = cell(length(y_ticks), 1);
for i = 1:length(y_ticks)
    % Las virgulillas (~) empujan el número a la izquierda
    y_labels{i} = sprintf('%d~~', y_ticks(i));
end
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'YTickLabel', y_labels);

% Etiquetas de los ejes alejadas con saltos de línea \n
xlabel(sprintf('X\n\n\n'));
ylabel(sprintf('\n\n\n Y'));
zlabel(sprintf('Z\n\n\n'));

view(-37.5, 60);

% Margen interno para que no se recorte nada al exportar a SVG
set(gca, 'Position', [0.15 0.15 0.7 0.8]);

hold off;

% 6. Guardar en formato SVG
print(NOMBRE_ARCHIVO_SVG, '-dsvg');
disp(['Gráfico final guardado en: ', fullfile(pwd, NOMBRE_ARCHIVO_SVG)]);
