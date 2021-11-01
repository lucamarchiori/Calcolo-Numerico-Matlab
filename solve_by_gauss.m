function [xs, L, R] = solve_by_gauss(A, b)
% GAUSS1 - Fattorizzazione di Gauss, versione 1
% L = Triangolare inferiore
% R = Triangolare superiore
% xs =
% A = Matrice coefficienti
% b = Vettore soluzioni

soglia_verso_zero = eps * norm(A,inf);
[mA, nA] = size(A);
[mb, nb] = size(b)

%Se le righe del vettore soluzione é diverso dal numero di righe di A
if(mb ~= ma)
    error('Fattorizzazione non effettuabile.');
end

Im = eye(mA); %Matrice identita
Ab = [A,b]; %Matrice completa
L = Im;
n_passi = min(mA-1,nA);
for k = 1 : n_passi-1
    % Controllo se elementi diagonali a zero
    if ( abs( A(k,k) ) <  soglia_verso_zero)
        error('Fattorizzazione non effettuabile.');
    else
        %mk = Vettore dei moltiplicatori. Vettore colonna con tante righe quante
        %quelle della matrice. Il primo elemento é zero, gli altri sono uguali
        % al rapporto tra l'elemento in posizione k e quello sulla diagonale
        mk=[zeros(k,1); A((k+1):end,k)/A(k,k)]; 
        %ekT = riga k-esima della matrice identità 
        ekT = Im(k,:);
        %Lk = matrice identità con il vettore dei moltiplicatori e tutto il
        %resto a zero
        Lk = Im - mk*ekT;
        %Ab = nuova matrice completa
        Ab = Lk * Ab;
        invLk = Im + mk*ekT;
        L = L*invLk;
    end
end
R = solve_tri_sup_by_row(Ab(:,1:nA)); %Passo solo la parte di A, escludo b
y = Ab(:,end);
xs = solve_tri_low_by_row(R,y);
end