% Parámetros base
Nel = 16;
p_values = [2, 3, 4, 5];

% Figura de altura 280 para acomodar textos arriba y abajo
fig = figure('Position', [100, 100, 1000, 280], 'Visible', 'off');

% Mapa de colores: Blanco para 0, Azul oscuro para 1
mi_mapa = [1, 1, 1; 0.1, 0.2, 0.4];
colormap(mi_mapa);

% Definición manual de márgenes (en unidades normalizadas 0 a 1)
margen_izq = 0.02;
margen_inf = 0.15;
espacio_entre = 0.02;
ancho_eje = (1 - 2*margen_izq - 3*espacio_entre) / 4;
% Reducimos de 0.72 a 0.70 para dar más espacio por arriba y evitar cortes
alto_eje = 0.70;

for k = 1:length(p_values)
    p = p_values(k);
    Nb = Nel + p;

    K = zeros(Nb, Nb);
    for i = 1:Nb
        j_min = max(1, i - p);
        j_max = min(Nb, i + p);
        K(i, j_min:j_max) = 1;
    end

    % --- Cálculo de la proporción de elementos no nulos ---
    elementos_no_nulos = sum(K(:));
    elementos_totales = Nb * Nb;
    porcentaje = (elementos_no_nulos / elementos_totales) * 100;

    % Calculamos la posición exacta del bloque k
    pos_x = margen_izq + (k-1) * (ancho_eje + espacio_entre);
    ax = axes('Position', [pos_x, margen_inf, ancho_eje, alto_eje]);

    imagesc(K);

    % Estética de la matriz
    axis square;
    axis ij;
    set(ax, 'XTick', [], 'YTick', []); % Eliminamos números de ejes

    % --- CORRECCIÓN DEL TÍTULO SUPERIOR ---
    % Usamos Position [0.5, 1.03] para acercarlo a la matriz y bajarlo
    title(sprintf('p = %d', p), 'FontSize', 12, 'FontWeight', 'bold', ...
          'Units', 'normalized', 'Position', [0.5, 1.03]);

    % --- TEXTO INFERIOR ---
    texto_inferior = sprintf('%d / %d (%.1f%%)', elementos_no_nulos, elementos_totales, porcentaje);
    text(0.5, -0.12, texto_inferior, 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', 11);

    % Cuadrícula ultra fina
    hold on;
    line_color = [0.85, 0.85, 0.85];
    for idx = 0.5 : 1 : Nb+0.5
        plot([0.5, Nb+0.5], [idx, idx], 'Color', line_color, 'LineWidth', 0.3);
        plot([idx, idx], [0.5, Nb+0.5], 'Color', line_color, 'LineWidth', 0.3);
    end
    hold off;

    xlim([0.5, Nb+0.5]);
    ylim([0.5, Nb+0.5]);
end

% Ajuste final: Fondo blanco puro
set(fig, 'Color', 'w');

% EXPORTACIÓN CRÍTICA:
print(fig, 'densidad_matriz_final.png', '-dpng', '-r300');

disp('Figura generada: densidad_matriz_final.png con título ajustado y proporciones visibles.');
