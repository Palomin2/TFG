% ==== Parámetros (ajusta aquí los que correspondan a tu fichero) ====
hX = 16;      % Ejemplo
pX = 2;      % Ejemplo
hY = 16;     % Ejemplo
pY = 2;      % Ejemplo

% ==== Construir ruta del archivo exactamente como en C++ ====
base_path = "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/";
filename = sprintf("surfaceCircle3d_hX=%d_pX=%d_hY=%d_p=%d data.txt", hX, pX, hY, pY);
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
figure;
surf(Xg, Yg, Zg, Vg);
shading interp;
colorbar;
xlabel('X');
ylabel('Y');
zlabel('Z');
title(sprintf('Superficie SolEvals (hX=%d, pX=%d, hY=%d, pY=%d)', hX, pX, hY, pY));
view(45,30);

