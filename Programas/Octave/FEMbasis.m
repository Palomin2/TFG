clear; clc; close all;

L = 1;
num_elementos = 6;
p = 1;

num_nodos = num_elementos * p + 1;
x_nodos = linspace(0, L, num_nodos);
colores = lines(num_nodos);

figure;
% --- AUMENTAMOS LOS MÁRGENES DE LA FIGURA ---
% [izquierda, abajo, ancho, alto] - Dejamos un 15% de margen por la izq y abajo
set(gca, 'Position', [0.15 0.15 0.8 0.8]);
% --------------------------------------------

hold on; grid on;
axis([-0.1 L+0.1 -0.2 1.2]);
box on; % Mantiene el contorno cerrado

% --- BLOQUE PARA DESPLAZAR LOS NÚMEROS DE LOS EJES ---
xticks_vals = linspace(0, L, 6);
yticks_vals = [0 0.2 0.4 0.6 0.8 1];

set(gca, 'XTick', xticks_vals);
set(gca, 'YTick', yticks_vals);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);

desplazamiento_y = -0.32; % Ajustado ligeramente para el nuevo margen
desplazamiento_x = -0.2;

for k = 1:length(xticks_vals)
    text(xticks_vals(k), desplazamiento_y, sprintf('%g', xticks_vals(k)), ...
        'HorizontalAlignment', 'center', 'Clipping', 'off');
end

for k = 1:length(yticks_vals)
    text(desplazamiento_x, yticks_vals(k), sprintf('%g', yticks_vals(k)), ...
        'HorizontalAlignment', 'right', 'Clipping', 'off');
end
% -----------------------------------------------------

for i = 1:num_nodos
    valores_nodales = zeros(1, num_nodos);
    valores_nodales(i) = 1;
    x_plot_total = [];
    y_plot_total = [];

    for e = 1:num_elementos
        idx_ini = (e-1)*p + 1;
        idx_fin = e*p + 1;
        indices_locales = idx_ini:idx_fin;
        x_elem = x_nodos(indices_locales);
        y_elem = valores_nodales(indices_locales);

        coeficientes = polyfit(x_elem, y_elem, p);
        xx_elem = linspace(x_elem(1), x_elem(end), 50);
        yy_elem = polyval(coeficientes, xx_elem);

        x_plot_total = [x_plot_total, xx_elem];
        y_plot_total = [y_plot_total, yy_elem];
    end

    plot(x_plot_total, y_plot_total, 'Color', colores(i,:), 'LineWidth', 2);
    plot(x_nodos(i), 1, 'o', 'MarkerFaceColor', colores(i,:), 'MarkerEdgeColor', 'k');
end

for e = 0:num_elementos
    x_borde = e * (L/num_elementos);
    line([x_borde x_borde], [-0.2 1.2], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
end

hold off;

nombre_archivo = sprintf('bases_fem_p%d.svg', p);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, nombre_archivo, '-dsvg');
