% ##########################################################################
% SCRIPT OCTAVE PARA TRAZAR EVALUACIÓN DE UN CAMPO SOBRE UNA CURVA NURBS
% Visualización: Dominio base (2D) + Elevación Z y Mapa de Color según Evaluación
% ##########################################################################

% --- 1. CONFIGURACIÓN DE RUTAS Y ARCHIVOS ---
RUTA_ARCHIVO_CURVA = 'C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Curvas\CurvaEjCircle_h=128_p=3.txt';
RUTA_ARCHIVO_EVAL  = 'C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConst_h=128_p=3_testCircle.txt';
NOMBRE_ARCHIVO_SVG = 'EvaluacionCurva_Colormap_Plot.svg';

disp('Cargando datos de geometría y evaluación...');

% --- 2. LECTURA Y PROCESAMIENTO DE DATOS ---
try
    D_curva = load(RUTA_ARCHIVO_CURVA);
    V_raw   = load(RUTA_ARCHIVO_EVAL);
catch
    error('ERROR: No se pudieron cargar los archivos. Comprueba los nombres y rutas.');
end

% Procesar geometría (Asegurar que sea N x 3)
[M, K] = size(D_curva);
if M == 3
    C = D_curva'; % Transponer de 3xN a Nx3
else
    C = D_curva;
end
X = C(:,1); Y = C(:,2); Z_base = C(:,3);
N_puntos = length(X);

% Procesar valores de evaluación (Asegurar que sea un vector columna N x 1)
V = V_raw(:);

% Verificación crítica de coherencia matemática
if length(V) ~= N_puntos
    error('DESAJUSTE DIMENSIONAL: La curva tiene %d puntos pero hay %d valores de evaluación.', N_puntos, length(V));
end

disp(['Datos validados: ', num2str(N_puntos), ' puntos cargados correctamente.']);

% --- 3. CONFIGURACIÓN DE LA FIGURA ---
figure('Position', [100, 100, 900, 700]);
hold on;
grid on;

% A. Trazar la curva dominio original (Circunferencia plana en Z=0)
plot3(X, Y, Z_base, 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.5);

% B. Trazar la curva evaluada (Z = Evaluación) con Mapa de Color
% Truco topológico en Octave para crear líneas con gradiente de color continuo:
% Se define una superficie de anchura cero mediante la duplicación de los vectores
surface([X, X]', [Y, Y]', [V, V]', [V, V]', ...
        'facecolor', 'none', ...
        'edgecolor', 'interp', ...
        'linewidth', 3);

% Configurar el mapa de color (Jet es el estándar para campos escalares, o 'parula')
colormap(jet);
cb = colorbar; % Guardamos el 'handle' de la barra para manipularla después
%title(cb, 'FontSize', 10, 'FontWeight', 'bold');

% --- 4. LÍMITES Y ESTÉTICA DE EJES (Calibración limpia) ---
% Extraer límites para ajustar la caja
min_X = floor(min(X)); max_X = ceil(max(X));
min_Y = floor(min(Y)); max_Y = ceil(max(Y));
min_V = floor(min(min(Z_base), min(V)));
max_V = ceil(max(max(Z_base), max(V)));

xlim([min_X - 0.2, max_X + 0.2]);
ylim([min_Y - 0.2, max_Y + 0.2]);
zlim([min_V - 0.2, max_V + 0.2]);

% Personalización de los 'Ticks' (marcas) como en tu script
set(gca, 'XTick', min_X:1:max_X);
set(gca, 'ZTick', -1:1:1);
set(gca, 'YTick', min_Y:1:max_Y);

% Estilo de marcas hacia afuera
set(gca, 'TickDir', 'out', 'TickLength', [0.1, 0.1]);

% Truco TeX para separación (eje Y)
y_labels = cell(length(min_Y:1:max_Y), 1);
for i = 1:length(min_Y:1:max_Y)
    y_labels{i} = sprintf('%d~~', min_Y - 1 + i);
end
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'YTickLabel', y_labels);

% Etiquetas de ejes espaciadas
xlabel(sprintf('X\n\n\n'));
ylabel(sprintf('\n\n\nY'));
zlabel(sprintf('Z\n\n\n'));

% Ángulo de vista óptimo para ver la base y la elevación
view(-37.5, 30);

% --- SEPARACIÓN FÍSICA DE LA BARRA DE COLOR ---
% Ajuste del marco interno de la gráfica [Izquierda, Abajo, Ancho, Alto]
% Reducimos el ancho a 0.65 (65%) para no chocar con la barra
set(gca, 'Position', [0.10 0.15 0.65 0.8]);

% Forzamos la posición de la barra de color separada a la derecha
% La ponemos en el 85% de la figura, dejando un 10% de espacio vacío entre ambas
set(cb, 'Position', [0.85 0.15 0.03 0.8]);

hold off;

% --- 5. EXPORTACIÓN ---
% Guardamos el SVG final
print(NOMBRE_ARCHIVO_SVG, '-dsvg');
disp(['Gráfica generada e imagen guardada como: ', NOMBRE_ARCHIVO_SVG]);
