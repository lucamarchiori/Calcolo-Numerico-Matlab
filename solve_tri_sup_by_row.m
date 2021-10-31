function [x] = solve_tri_sup_by_row(R,b)
%solve_tri_sup_by_row - Soluzione di sistema triangolare superiore (per righe)
%Algoritmo di sostituzione all'indietro.

%Si ricava x(n) dall'ultima equazione, si sostituisce alla penultima e si
%ricava x-1 e si procede iterativamente fino alla soluzione completa.

n = length(b);
x=b;
x(n) = x(n)/R(n,n);
for i=n-1:-1:1
    x(i) = x(i)-R(i,i+1:n)*x(i+1:n);
    x(i) = x(i)/R(i,i);
end
end