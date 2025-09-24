function [x,w] = gauss_legendre(n, a, b)
  % Devuelve n puntos y pesos de Gauss–Legendre en el intervalo [a,b]
  % usando la función legpts de Chebfun si existe, o eigenmétodo
  beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
  T = diag(beta,1) + diag(beta,-1);
  [V,D] = eig(T);
  x = diag(D);          % puntos en (-1,1)
  [x, i] = sort(x);
  w = 2*V(1,i).^2;      % pesos en (-1,1)
  % reescala a [a,b]
  x = 0.5*((b-a)*x + (b+a));
  w = 0.5*(b-a)*w;
end
