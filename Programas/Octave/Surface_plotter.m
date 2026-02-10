% ==== Parámetros (ajusta aquí los que correspondan a tu fichero) ====
hX = 256;      % Ejemplo
pX = 2;      % Ejemplo
hY = 256;     % Ejemplo
pY = 2;      % Ejemplo

% ==== Construir ruta del archivo exactamente como en C++ ====
base_path = "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/";
filename = sprintf("surfaceTrivial3d_hX=%d_pX=%d_hY=%d_p=%d data.txt", hX, pX, hY, pY);
fullpath = strcat(base_path, filename);
fprintf("Leyendo datos desde: %s\n", fullpath);
% ==== Cargar datos (X Y Z Valor) ====
data = load(fullpath);
X = data(:,1);
Y = data(:,2);
Z = data(:,4);
V = data(:,4);
% ==== Calcular dimensiones teóricas ====
nx = hX * (pX + 1);
ny = hY * (pY + 1);
% ==== Reorganizar en malla ====
Xg = reshape(X, nx, ny);
Yg = reshape(Y, nx, ny);
Zg = reshape(Z, nx, ny);
Vg = reshape(V, nx, ny);
% ==== Dibujar ====
hFig = figure;
surf(Xg, Yg, Zg, Vg);
shading interp;
colorbar;
xlabel('X');
ylabel('Y');
zlabel('Z');
%title(sprintf('Superficie SolEvals (hX=%d, pX=%d, hY=%d, pY=%d)', hX, pX, hY, pY));
view(15,20);

%Axis tight quita el espacio blanco extra
%axis tight;

% ... (Tu código de surf, view y axis manual) ...

% ==== TRUCO PARA QUITAR MÁRGENES BLANCOS ====

% 1. Poner el fondo blanco (por si acaso sale gris)
set(hFig, 'Color', 'w');
set(gca, 'Color', 'w');

% 2. Expandir los ejes para que ocupen TODO el espacio de la figura
%    'OuterPosition' incluye los números y etiquetas.
%    Al ponerlo a [0 0 1 1], los pegamos al borde de la imagen.
set(gca, 'OuterPosition', [0, 0, 1, 1]);

% 3. Reducir el margen interno "de seguridad" que pone Octave
%    (TightInset es el espacio mínimo necesario para las letras)
ti = get(gca, 'TightInset');
set(gca, 'Position', [ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

% 4. Configurar el modo de papel para que respete lo que ves en pantalla
set(hFig, 'PaperPositionMode', 'auto');

% ==== GUARDAR EN PNG DE ALTA CALIDAD (600 DPI) ====
% Usamos -r600 para que al ponerlo en LaTeX al 50% de tamaño
% tenga una densidad de píxeles brutal (como una Retina Display).
output_filename = sprintf("surface_hX=%d_CROP_HQ1.png", hX);
output_fullpath = strcat(base_path, output_filename);

fprintf("Guardando PNG recortado HQ en: %s\n", output_fullpath);

% drawnow asegura que los cambios de márgenes se apliquen antes de imprimir
drawnow;
print(hFig, output_fullpath, "-dpng", "-r600");
