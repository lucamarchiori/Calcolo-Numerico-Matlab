function [x] = solve_tri_low_by_row(L, b)
% solve_tri_low_by_row - Soluzione di sistemi triangolari inferiori (per righe)
% Algoritmo di eliminazione in avanti

% Per i sistemi triangolari inferiori si usa l algoritmo di eliminazione in avanti,
% ove si ricava x1 dalla prima equazione, si sostituisce nella seconda e si ricava
% x2 e così via (algoritmo per righe).
n = length(b);
x = b;
x(1) = x(1)/L(1,1);
for i = 2:n
    x(i) = x(i) - L(i, 1:i-1)*x(1:i-1);
    x(i) = x(i)/L(i,i);
end
end