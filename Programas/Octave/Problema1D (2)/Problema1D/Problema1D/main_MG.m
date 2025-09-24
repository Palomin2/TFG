%% ALGORITMO DE MULTIGRID

tol = 10^(-30);
smth=1;
w = 1;
cyc_ind = 2;  % 1 for V-cycle, and 2 for W-cycle
R = cell(max_level,1);
A = cell(max_level,1);
P = cell(max_level,1);
%----------------
% finest level
%----------------
A{1} = Af;
h2 = h;
Nc = N;
for level = 2:max_level
    Nc = Nc/2;
    h2 = 2*h2; 
    R{level} = restriction(Nc);
    P{level} = 4*R{level}';
    A{level} = R{level}*A{level-1}*P{level};
end
res = 1.0;
res_old = 1.0;
res_in_p = max(abs(bp-Af*p)); 
iter = 1;
%while res>tol*res_in_p   
for iter = 1:30
      p = MG(A,R,P,p,bp,h,nu1,nu2,N,w,cyc_ind,smth,max_level,1);      
      rp = bp - Af*p;   
      res = max(abs(rp));
      RR = reshape(rp,N,N);
%       iter
       res/res_old
      factor(iter) = res/res_old;      
      res_old = res;
end
rho_m = mean(factor(20:30))