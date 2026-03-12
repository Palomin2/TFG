function plot_bernstein_bases(p)
    m = 100;
    u = linspace(0, 1, m+1);

    figure;

    % --- AUMENTAMOS LOS MÁRGENES DE LA FIGURA ---
    % Dejamos un 15% de margen por la izquierda y abajo para que entren los textos
    set(gca, 'Position', [0.15 0.15 0.8 0.8]);
    % --------------------------------------------

    hold on; grid on;

    % Damos aire a los límites para que las curvas no toquen los bordes
    axis([-0.1 1.1 -0.2 1.2]);
    box on; % Mantiene el contorno cerrado

    % --- BLOQUE PARA DESPLAZAR LOS NÚMEROS DE LOS EJES ---
    xticks_vals = 0:0.2:1; % Marcas cada 0.2 en X
    yticks_vals = 0:0.2:1; % Marcas cada 0.2 en Y

    set(gca, 'XTick', xticks_vals);
    set(gca, 'YTick', yticks_vals);

    % Vaciamos las etiquetas automáticas de Octave
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);

    % Coordenadas manuales fuera de la caja para los números
    desplazamiento_y = -0.32;
    desplazamiento_x = -0.2;

    % Dibujamos los números del eje X a mano
    for k = 1:length(xticks_vals)
        text(xticks_vals(k), desplazamiento_y, sprintf('%g', xticks_vals(k)), ...
            'HorizontalAlignment', 'center', 'Clipping', 'off');
    end

    % Dibujamos los números del eje Y a mano
    for k = 1:length(yticks_vals)
        text(desplazamiento_x, yticks_vals(k), sprintf('%g', yticks_vals(k)), ...
            'HorizontalAlignment', 'right', 'Clipping', 'off');
    end
    % -----------------------------------------------------

    colors = hsv(p + 1);

    for i = 0:p
        binom_coeff = nchoosek(p, i);
        poly_term = u.^i .* (1 - u).^(p - i);
        B_ip = binom_coeff * poly_term;

        plot(u, B_ip, 'LineWidth', 2, 'Color', colors(i + 1, :));
    end

    hold off;

    % --- EXPORTACIÓN ---
    filename = sprintf('bernstein_bases_p%d.svg', p);
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, filename, '-dsvg');

    disp(['Figura guardada como: ', filename]);
end
