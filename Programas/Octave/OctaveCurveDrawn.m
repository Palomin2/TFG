%% ============================================================
%% Script: visualizar_curva_coloreada.m
%% Descripción:
%%   Lee una curva en R^3 desde archivo y un vector de valores
%%   escalares asociados, y genera una visualización 3D coloreada.
%% ============================================================

clear; close all; clc;

%% --- Parámetros configurables ---
h_val = 256;
p_val = 2;
test_name = 'testCircle';

%% --- Archivos de entrada ---
% Ajusta las rutas según tus necesidades:
file_curve = sprintf('C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/CurvaEjCircle_h=%d_p=%d.txt', h_val, p_val);
file_values = sprintf('C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/EjSinNonConst_h=%d_p=%d_Analytic_%s.txt', h_val, p_val, test_name);

%% --- Leer datos ---
% La curva tiene tres filas: X, Y, Z
curve_data = dlmread(file_curve);
x = curve_data(1, :);
y = curve_data(2, :);
z = curve_data(3, :);

% Vector de valores asociados
values = dlmread(file_values);
values = values(:)'; % asegurar forma de fila

n_points = length(values);
if length(x) != n_points
    error('La curva y el vector de valores no tienen el mismo número de puntos.');
endif

%% --- Figura 1: Curva coloreada ---
h1 = figure('Position', [100 100 1000 800]);
hold on; grid on; axis equal;

% Graficar con color dependiente del valor escalar
colormap jet; % puedes usar 'parula', 'turbo', 'hot', etc.
scatter3(x, y, values, 40, values, 'filled'); % puntos coloreados
plot3(x, y, values, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0); % línea gris

colorbar;
%title(sprintf('Curva en R^3 coloreada (h=%d, p=%d)', h_val, p_val), 'FontSize', 20);
xlabel('x'); ylabel('y'); zlabel('z');
view(3); % vista 3D

% Guardar figura
print(h1, sprintf('CurvaColor_h%d_p%d_Z.svg', h_val, p_val), '-dsvg');
print(h1, sprintf('CurvaColor_h%d_p%d_Z.png', h_val, p_val), '-dpng', '-r300');
hold off;

fprintf('Figura 3D generada y guardada correctamente.\n');
