% Script di prova degli algoritmi di fattorizzazione LDU di Gauss

A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1];
b = [4; 1; -3; 4];

[soluzione, L, R] = solve_by_gauss(A,b)
detA = prod(diag(R))