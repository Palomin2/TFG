% ##########################################################################
% SCRIPT OCTAVE PARA TRAZAR EVALUACIÓN DE UN CAMPO SOBRE UNA CURVA NURBS
% Visualización: Curva en su Z original (Plana) + Mapa de Color según Evaluación
% ##########################################################################

% --- 1. CONFIGURACIÓN DE RUTAS Y ARCHIVOS ---
RUTA_ARCHIVO_CURVA = 'C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\Curvas\CurvaEjCircle_h=128_p=3.txt';
RUTA_ARCHIVO_EVAL  = 'C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\SolEvals\EjSinNonConst_h=128_p=3_testCircle.txt';
NOMBRE_ARCHIVO_SVG = 'EvaluacionCurva_Plana_Colormap.svg';

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

% Procesar valores de evaluación
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

% Trazar la curva en su posición Z normal, pero coloreada usando V
% Argumentos de surface: X, Y, Z, Color
surface([X, X]', [Y, Y]', [Z_base, Z_base]', [V, V]', ...
        'facecolor', 'none', ...
        'edgecolor', 'interp', ...
        'linewidth', 3);

% Configurar el mapa de color
colormap(jet);
cb = colorbar;

% --- 4. LÍMITES Y ESTÉTICA DE EJES ---
% Extraer límites X e Y
min_X = floor(min(X)); max_X = ceil(max(X));
min_Y = floor(min(Y)); max_Y = ceil(max(Y));

% Manejo seguro de Z para curvas planas (como Z=0)
min_Z = floor(min(Z_base));
max_Z = ceil(max(Z_base));
if min_Z == max_Z
    min_Z = min_Z - 1;
    max_Z = max_Z + 1;
end

xlim([min_X - 0.2, max_X + 0.2]);
ylim([min_Y - 0.2, max_Y + 0.2]);
zlim([min_Z, max_Z]);

% Personalización de los 'Ticks'
set(gca, 'XTick', min_X:1:max_X);
set(gca, 'YTick', min_Y:1:max_Y);
set(gca, 'ZTick', min_Z:1:max_Z);

% Estilo de marcas hacia afuera
set(gca, 'TickDir', 'out', 'TickLength', [0.1, 0.1]);

% Truco TeX para separación (eje Y)
y_labels = cell(length(min_Y:1:max_Y), 1);
for i = 1:length(min_Y:1:max_Y)
    y_labels{i} = sprintf('%d~~', min_Y - 1 + i);
end
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'YTickLabel', y_labels);

% Etiquetas de ejes
xlabel(sprintf('X\n\n\n'));
ylabel(sprintf('\n\n\nY'));
zlabel(sprintf('Z\n\n\n'));

% Ángulo de vista (Mantenido en 3D para notar que está en Z=0)
% NOTA: Si prefieres verla totalmente desde arriba en 2D puro, cambia a view(0, 90);
view(-37.5, 30);

% --- SEPARACIÓN FÍSICA DE LA BARRA DE COLOR ---
% Caja de los ejes 3D ocupando hasta el 65% del ancho
set(gca, 'Position', [0.10 0.15 0.65 0.8]);

% Barra de color empujada al 85% de la ventana
set(cb, 'Position', [0.85 0.15 0.03 0.8]);

hold off;

% --- 5. EXPORTACIÓN ---
print(NOMBRE_ARCHIVO_SVG, '-dsvg');
disp(['Gráfica generada e imagen guardada como: ', NOMBRE_ARCHIVO_SVG]);
