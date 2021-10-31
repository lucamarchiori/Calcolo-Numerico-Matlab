function [x] = solve_tri_sup_by_col(R,b)
%solve_tri_sup_by_col - Soluzione di sistema triangolare superiore (per colonne)
%Algoritmo di sostituzione all'indietro.

% Dopo aver ricavato x(n) dall' ultima equazione, sostituire la soluzione in tutte le 
% le equazioni precedenti ottenendo un sistema triangolare superiore di dimensione n-1;
% da questo si ricava l'ultima incognita x(n-1) dall'ultima equazione e si 
% sostituisce nelle precedenti e così via.

n = length(b);
x=b;
x(n) = x(n)/R(n,n);
for j=n-1:-1:1
    x(1:j) = x(1:j)-R(1:j,j+1)*x(j+1);
    x(j) = x(j)/R(j,j);
end
end