clear; clc; close all;

% --- 1. PARÁMETROS CONFIGURABLES ---
L = 1;
num_elementos = 6;
p = 1;

% --- 2. GENERACIÓN DE LA MALLA ---
num_nodos = num_elementos * p + 1;
x_nodos = linspace(0, L, num_nodos);
colores = lines(num_nodos);

% Configuración de la figura
hold on; grid on;
axis([-0.1 L+0.1 -0.2 1.2]);

% --- 3. BUCLE DE FUNCIONES BASE ---
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
    % --- DIBUJO ---
    plot(x_plot_total, y_plot_total, 'Color', colores(i,:), 'LineWidth', 2);
    %try
        %fill(x_plot_total, y_plot_total, colores(i,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    %catch
    %end
    plot(x_nodos(i), 1, 'o', 'MarkerFaceColor', colores(i,:), 'MarkerEdgeColor', 'k');
    if p == 1
        offset = 0.05;
    else
        offset = 0.1 + 0.05*mod(i,2);
    end
    text(x_nodos(i), 1 + offset, sprintf('\\phi_{%d}', i), ...
         'Color', colores(i,:), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
for e = 0:num_elementos
    x_borde = e * (L/num_elementos);
    line([x_borde x_borde], [-0.2 1.2], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
end

hold off;

nombre_archivo = sprintf('bases_fem_p%d.svg', p);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, nombre_archivo, '-dsvg');
