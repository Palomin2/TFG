% --- CARGA DE DATOS ---
ruta_bases = 'C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\RatBasisFuncts.txt';
ruta_derivadas = 'C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\DataFiles\Nurbs\RatBasisFuncts_Derivates.txt';

try
    data1 = dlmread(ruta_bases);
    data2 = dlmread(ruta_derivadas);
catch
    error('ERROR: No se pudieron cargar los archivos de datos. Revisa las rutas.');
end

[n1, n2] = size(data1);
n1=n1-1;
% Vector paramétrico normalizado [0, 1]
u = linspace(0, 1, n2);

% Paleta HSV para colores vivos
colores = hsv(n1) * 0.93;

% ##########################################################################
% FIGURA 1: FUNCIONES BASE NURBS
% ##########################################################################
figure(1);
clf;
% Aumentamos el margen (de 0.15 a 0.18) para que quepan los números más alejados
set(gca, 'Position', [0.15 0.15 0.75 0.75]);
hold on; grid on; box on;

axis([-0.1 1.1 -0.2 1.2]);

xticks_vals = 0:0.2:1;
yticks_vals = 0:0.2:1;
set(gca, 'XTick', xticks_vals, 'YTick', yticks_vals);
set(gca, 'XTickLabel', [], 'YTickLabel', []);

% --- AJUSTES DE POSICIÓN FIGURA 1 ---
desplazamiento_y = -0.35; % Ligeramente más abajo
desplazamiento_x = -0.22; % MÁS A LA IZQUIERDA (antes era -0.15)
% ------------------------------------

for k = 1:length(xticks_vals)
    text(xticks_vals(k), desplazamiento_y, sprintf('%g', xticks_vals(k)), ...
        'HorizontalAlignment', 'center', 'Clipping', 'off', 'FontSize', 12);
end
for k = 1:length(yticks_vals)
    text(desplazamiento_x, yticks_vals(k), sprintf('%g', yticks_vals(k)), ...
        'HorizontalAlignment', 'right', 'Clipping', 'off', 'FontSize', 12);
end

for i = 1:n1
    plot(u, data1(i,:), 'LineWidth', 2.5, 'Color', colores(i,:));
end
hold off;

print('funciones_base_NURBS_3.svg', '-dsvg');
disp('Figura 1 guardada como: funciones_base_NURBS.svg');


% ##########################################################################
% FIGURA 2: DERIVADAS DE LAS FUNCIONES BASE
% ##########################################################################
figure(2);
clf;
% Aumentamos el margen aquí también
set(gca, 'Position', [0.18 0.18 0.75 0.75]);
hold on; grid on; box on;

min_der = min(data2(:));
max_der = max(data2(:));
rango_y = max_der - min_der;
margen_y = rango_y * 0.15;

lim_inf = min_der - margen_y;
lim_sup = max_der + margen_y;
axis([-0.1 1.1 lim_inf lim_sup]);

yticks_der = linspace(min_der, max_der, 5);
set(gca, 'XTick', xticks_vals, 'YTick', yticks_der);
set(gca, 'XTickLabel', [], 'YTickLabel', []);

% --- AJUSTES DE POSICIÓN FIGURA 2 ---
% Bajamos el eje X multiplicando por 0.14 (antes era 0.08)
desp_y_der = lim_inf - (rango_y * 0.14);
% El eje Y usa el mismo desplazamiento_x = -0.22 de la Figura 1
% ------------------------------------

for k = 1:length(xticks_vals)
    text(xticks_vals(k), desp_y_der, sprintf('%g', xticks_vals(k)), ...
        'HorizontalAlignment', 'center', 'Clipping', 'off', 'FontSize', 12);
end
for k = 1:length(yticks_der)
    text(desplazamiento_x, yticks_der(k), sprintf('%.1f', yticks_der(k)), ...
        'HorizontalAlignment', 'right', 'Clipping', 'off', 'FontSize', 12);
end

for i = 1:n1
    plot(u, data2(i,:), 'LineWidth', 2.5, 'Color', colores(i,:));
end
hold off;

print('Ders_funciones_base_NURBS.svg', '-dsvg');
disp('Figura 2 guardada como: Ders_funciones_base_NURBS.svg');
