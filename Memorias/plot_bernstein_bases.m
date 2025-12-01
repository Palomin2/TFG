function plot_bernstein_bases(p)

    m = 100;
    u = linspace(0, 1, m+1);

    figure;
    hold on;
    grid on;

    colors = hsv(p + 1);

    for i = 0:p


        binom_coeff = nchoosek(p, i);
        poly_term = u.^i .* (1 - u).^(p - i);
        B_ip = binom_coeff * poly_term;
        plot(u, B_ip, 'LineWidth', 2, 'Color', colors(i + 1, :), ...
             'DisplayName', ['B_{', num2str(i), ',', num2str(p), '}(u)']);
    end

    xlim([0, 1]);
    ylim([0, 1.1]);
    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXPORTACIÓN EN FORMATO SVG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    filename = ['bernstein_bases_p', num2str(p), '.svg'];
    print(filename, '-dsvg');

    disp(['Figura guardada como: ', filename]);
end
