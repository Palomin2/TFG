% --- Test de integrales ---
[a,b] = deal(0,1);
n = 10;  % número de puntos de Gauss
[xq, wq] = gauss_legendre(n,a,b);

% Funciones de prueba
f1 = @(x) 1;
f2 = @(x) x;
f3 = @(x) x.^2;
f4 = @(x) x.^3;
f5 = @(x) x.^4;

% Integrales aproximadas
I1 = sum(wq .* f1(xq))
I2 = sum(wq .* f2(xq))
I3 = sum(wq .* f3(xq))
I4 = sum(wq .* f4(xq))
I5 = sum(wq .* f5(xq))

% Integrales exactas
I1_exact = 1;
I2_exact = 0.5;
I3_exact = 1/3;
I4_exact = 1/4;
I5_exact = 1/5;

% Mostrar errores
fprintf('Errores: \n');
fprintf('f1: %g\n', abs(I1-I1_exact));
fprintf('f2: %g\n', abs(I2-I2_exact));
fprintf('f3: %g\n', abs(I3-I3_exact));
fprintf('f4: %g\n', abs(I4-I4_exact));
fprintf('f5: %g\n', abs(I5-I5_exact));

