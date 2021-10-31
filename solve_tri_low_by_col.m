function [x] = solve_tri_low_by_col(L, b)
% solve_tri_low_by_col - Soluzione di sistemi triangolari inferiori (per colonne)
% Algoritmo di eliminazione in avanti

% Dopo aver ricavato x1 dalla prima equazione, sostituire in tutte le successive 
% ottenendo un sistema triangolare inferiore di dimensione n-1; da questo si 
% ricava la prima incognita x2 dalla prima equazione e si sostituisce nelle 
% successive equazioni e così via (algoritmo di eliminazione in avanti per colonne).
n = length(b);
x = b;
x(1) = x(1)/L(1,1);
for j = 2:n
    x(j:n) = x(j:n) - L(j:n, j-1)*x(j-1);
    x(j) = x(j) / L(j,j);
end
end