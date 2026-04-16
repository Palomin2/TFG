% Script para graficar y exportar las bases NURBS directamente a PNG
clear all; clc; clf;

% 1. Parametros
p = 1; U = [0, 0, 1, 1];
q = 2; V = [0, 0, 0, 1, 1, 1];

Px = [1, 1, 0;
      3, 3, 0];
Py = [0, 1, 1;
      0, 3, 3];
W  = [1, sqrt(2)/2, 1;
      1, sqrt(2)/2, 1];

num_pts = 40;
u_vec = linspace(0, 1, num_pts);
v_vec = linspace(0, 1, num_pts);

% 2. Funcion de Cox-de Boor
function N = cox_de_boor(i, d, t, knots)
    if d == 0
        if (t >= knots(i) && t < knots(i+1)) || (t == knots(end) && t == knots(i+1) && knots(i) ~= knots(i+1))
            N = 1.0;
        else
            N = 0.0;
        end
    else
        den1 = knots(i+d) - knots(i);
        den2 = knots(i+d+1) - knots(i+1);
        term1 = 0.0; term2 = 0.0;
        if den1 > 0, term1 = ((t - knots(i)) / den1) * cox_de_boor(i, d-1, t, knots); end
        if den2 > 0, term2 = ((knots(i+d+1) - t) / den2) * cox_de_boor(i+1, d-1, t, knots); end
        N = term1 + term2;
    end
end

% 3. Inicializacion de matrices
X_mat = zeros(num_pts, num_pts);
Y_mat = zeros(num_pts, num_pts);
Z_mat = cell(2, 3);
for i=1:2
    for j=1:3
        Z_mat{i,j} = zeros(num_pts, num_pts);
    end
end

% 4. Evaluacion
for j_idx = 1:num_pts
    v = v_vec(j_idx);
    for i_idx = 1:num_pts
        u = u_vec(i_idx);

        x_val = 0.0; y_val = 0.0; w_val = 0.0;
        Nu = zeros(2, 1); Mv = zeros(3, 1);

        for ii = 1:2, Nu(ii) = cox_de_boor(ii, p, u, U); end
        for jj = 1:3, Mv(jj) = cox_de_boor(jj, q, v, V); end

        for ii = 1:2
            for jj = 1:3
                basis = Nu(ii) * Mv(jj) * W(ii,jj);
                w_val = w_val + basis;
                x_val = x_val + basis * Px(ii,jj);
                y_val = y_val + basis * Py(ii,jj);
            end
        end

        X_mat(i_idx, j_idx) = x_val / w_val;
        Y_mat(i_idx, j_idx) = y_val / w_val;

        % Guardar el valor Z para cada una de las 6 bases
        for ii = 1:2
            for jj = 1:3
                Z_mat{ii,jj}(i_idx, j_idx) = (Nu(ii) * Mv(jj) * W(ii,jj)) / w_val;
            end
        end
    end
end

% 5. Graficar
figure(1, 'Position', [100, 100, 800, 600]); % Ventana mas grande para mejor resolucion
hold on; grid on;

% Paleta de colores similar a la que ibamos a usar en LaTeX (Azul, Rojo, Verde, Naranja, Cyan, Purpura)
colors = [0 0.4470 0.7410;  0.8500 0.3250 0.0980;  0.4660 0.6740 0.1880;
          0.9290 0.6940 0.1250;  0.3010 0.7450 0.9330;  0.4940 0.1840 0.5560];
c_idx = 1;

for ii = 1:2
    for jj = 1:3
        surf(X_mat, Y_mat, Z_mat{ii,jj}, ...
             'FaceColor', colors(c_idx, :), ...
             'EdgeColor', [0.2 0.2 0.2], ...
             'EdgeAlpha', 0.4, ...
             'FaceAlpha', 0.6);
        c_idx = c_idx + 1;
    end
end

view(40, 35);
xlabel('x'); ylabel('y'); zlabel('\phi_{i,j}(x,y)', 'FontSize', 14');
axis tight;
set(gca, 'FontSize', 12); % Tamaño de fuente legible

% 6. Guardar imagen
print('nurbs_bases_completas.png', '-dpng', '-r300');
disp('Imagen nurbs_bases_completas.png generada con exito en el directorio actual.');
