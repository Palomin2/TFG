RUTA_ARCHIVO_CURVA = PATH + 'EjCurva3d.txt';
NOMBRE_ARCHIVO_SVG = 'CurvaNurbs_3D_Plot.svg';

PC = [
    10, 0, 0;
    15, 10, 5;
    5, 15, 20;
    0, 5, 10;
    10, 0, 0
];

try
    D = load(RUTA_ARCHIVO_CURVA);
catch
    error('ERROR: No se pudo cargar el archivo. Comprueba que la ruta sea correcta: %s', RUTA_ARCHIVO_CURVA);
end

[M, N] = size(D);

C = D';

figure;
hold on;
grid on;


plot3(C(:,1), C(:,2), C(:,3), 'b-', 'LineWidth', 2);
plot3(PC(:,1), PC(:,2), PC(:,3), 'r--', 'LineWidth', 1, 'DisplayName', 'Polígono de Control');
plot3(PC(:,1), PC(:,2), PC(:,3), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', 'Puntos de Control');

min_coords = min([C; PC]);
max_coords = max([C; PC]);
margen = 2;
axis([min_coords(1)-margen max_coords(1)+margen ...
      min_coords(2)-margen max_coords(2)+margen ...
      min_coords(3)-margen max_coords(3)+margen]);

view(3);
set(gca, 'Box', 'on');
print(NOMBRE_ARCHIVO_SVG, '-dsvg');
disp(['Gráfico final guardado en: ', fullfile(pwd, NOMBRE_ARCHIVO_SVG)]);
