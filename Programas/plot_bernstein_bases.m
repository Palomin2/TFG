function plot_rational_bernstein_bases(p)
    m = 100;
    u = linspace(0, 1, m+1);

    % --- VECTOR DE PESOS ---
    % Debe tener longitud p+1.
    % Por defecto lo inicializamos a 1 (lo que devuelve la base estándar).
    w = [1, 10, 7, 1];

    % EJEMPLO: Si p>=2, alteramos un peso interno para romper la simetría
    % y observar el efecto de "atracción" racional en la gráfica.
    if p >= 2
        w(2) = 5; % Modifica este vector libremente según tus pruebas
    end
    % -----------------------

    figure;

    % --- AUMENTAMOS LOS MÁRGENES DE LA FIGURA ---
    set(gca, 'Position', [0.15 0.15 0.8 0.8]);
    % --------------------------------------------

    hold on; grid on;
    axis([-0.1 1.1 -0.2 1.2]);
    box on;

    % --- BLOQUE PARA DESPLAZAR LOS NÚMEROS DE LOS EJES ---
    xticks_vals = 0:0.2:1;
    yticks_vals = 0:0.2:1;

    set(gca, 'XTick', xticks_vals);
    set(gca, 'YTick', yticks_vals);

    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);

    desplazamiento_y = -0.32;
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

    colors = hsv(p + 1);

    % FASE 1: Calcular B_ip estándar y ensamblar el denominador W(u)
    B = zeros(p + 1, length(u));
    W = zeros(1, length(u));

    for i = 0:p
        binom_coeff = nchoosek(p, i);
        poly_term = u.^i .* (1 - u).^(p - i);
        B(i + 1, :) = binom_coeff * poly_term;

        % Sumatorio para la función de peso global W(u)
        W = W + B(i + 1, :) * w(i + 1);
    end

    % FASE 2: Calcular R_ip racional y dibujar
    for i = 0:p
        R_ip = (B(i + 1, :) * w(i + 1)) ./ W;
        plot(u, R_ip, 'LineWidth', 2, 'Color', colors(i + 1, :));
    end

    hold off;

    % --- EXPORTACIÓN ---
    filename = sprintf('rational_bernstein_bases_p%d_weighted.svg', p);
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, filename, '-dsvg');

    disp(['Figura guardada como: ', filename]);
end
