
m = 100;
u = 0:1/m:1;

P = [0, 1, 2, 3;
     0, 5, 0, 4.5;
     0, 1, 2, 4];

W = [1.0, 10.0, 7.0, 1.0];

[aux, n_puntos] = size(P);

sol = zeros(aux, m+1);

for k = 1:m+1

    Bold = zeros(n_puntos, 1);
    Bold(1) = 1;

    for j = 2:n_puntos
        Bnew = zeros(n_puntos, 1);
        for i = 1:j-1
            % B(i,j) = (1-u)*B(i,j-1) + u*B(i-1,j-1)
            Bnew(i) = Bnew(i) + (1 - u(k)) * Bold(i);
            Bnew(i+1) = Bnew(i+1) + u(k) * Bold(i);
        endfor
        Bold = Bnew;
    endfor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RACIONALIZACIÓN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Numerator = zeros(aux, 1);
    Denominator = 0;

    for j = 1:n_puntos

        WeightedBasis = Bold(j) * W(j);
        Numerator = Numerator + WeightedBasis * P(:, j);
        Denominator = Denominator + WeightedBasis;
    endfor

    sol(:,k) = Numerator / Denominator;

endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTEO Y AJUSTE DE EJES DINÁMICOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_coords = min([P, sol], [], 2);
max_coords = max([P, sol], [], 2);

padding = 0.2;

min_X = min_coords(1) - padding;
max_X = max_coords(1) + padding;
min_Y = min_coords(2) - padding;
max_Y = max_coords(2) + padding;
min_Z = min_coords(3) - padding;
max_Z = max_coords(3) + padding;

figure
hold on

% Plotear el polígono de control (rojo discontinuo con círculos)
plot3(P(1,:), P(2,:), P(3,:), 'r--o', 'LineWidth', 1, 'MarkerSize', 5);

% Plotear la curva racional (azul sólido)
plot3(sol(1,:), sol(2,:), sol(3,:), 'b', 'LineWidth', 2);

% ** Establecer límites dinámicos **
xlim([min_X, max_X]);
ylim([min_Y, max_Y]);
zlim([min_Z, max_Z]);

xlabel('X');
ylabel('Y');
zlabel('Z');
grid on
view(3);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORTACIÓN DE LA FIGURA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print('curva_bezier_racional.svg', '-dsvg');
