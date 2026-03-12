%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINICIÓN DE DATOS Y CÁLCULO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 100;
u = 0:1/m:1;

P = [0, 1, 2, 3;
     0, 5, 0, 4.5;
     0, 1, 2, 4];

W = [1.0, 1.0, 1.0, 1.0];

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

% Redondear a enteros para tener ejes limpios
min_X = floor(min_coords(1));
max_X = ceil(max_coords(1));
min_Y = floor(min_coords(2));
max_Y = ceil(max_coords(2));
min_Z = floor(min_coords(3));
max_Z = ceil(max_coords(3));

% Crear figura con un tamaño específico
figure('Position', [100, 100, 800, 600])
hold on

% Plotear el polígono de control
plot3(P(1,:), P(2,:), P(3,:), 'r--o', 'LineWidth', 1, 'MarkerSize', 5);

% Plotear la curva racional
plot3(sol(1,:), sol(2,:), sol(3,:), 'b', 'LineWidth', 2);

% Establecer límites cerrados a números enteros
xlim([min_X, max_X]);
ylim([min_Y, max_Y]);
zlim([min_Z, max_Z]);

% Forzar las marcas
set(gca, 'XTick', min_X:1:max_X);
y_ticks = min_Y:1:max_Y; % Mantenemos el salto de 2 en 2
set(gca, 'YTick', y_ticks);
set(gca, 'ZTick', min_Z:1:max_Z);

% Volvemos a un TickLength prudente que no estorbe
set(gca, 'TickDir', 'out', 'TickLength', [0.1, 0.1]);

% ** TRUCO DEFINITIVO CON TeX: Espacios insecables a la derecha **
y_labels = cell(length(y_ticks), 1);
for i = 1:length(y_ticks)
    % Las virgulillas (~) actúan como espacios de hormigón que el SVG no puede borrar.
    % Cuantas más pongas, más a la izquierda se moverá el número.
    y_labels{i} = sprintf('%d~~', y_ticks(i));
end

% Asegurarnos de que Octave interpreta las ~ como espacios TeX
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'YTickLabel', y_labels);

% Tus etiquetas intactas
xlabel(sprintf('X\n\n\n'));
ylabel(sprintf('\n\n\n Y'));
zlabel(sprintf('Z\n\n\n'));

grid on
view(3);

% Margen interno
set(gca, 'Position', [0.15 0.15 0.7 0.8]);

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORTACIÓN DE LA FIGURA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print('curva_bezier.svg', '-dsvg');
