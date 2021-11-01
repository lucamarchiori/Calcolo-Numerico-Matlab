function [xs, L, R] = solve_by_gauss(A, b)
% GAUSS1 - Fattorizzazione di Gauss, versione 1
% L = Triangolare inferiore
% R = Triangolare superiore
% xs =
% A = Matrice coefficienti
% b = Vettore soluzioni

soglia_verso_zero = eps * norm(A,inf);
[mA, nA] = size(A);
[mb, nb] = size(b);

Im = eye(mA); %Matrice identita
Ab = [A,b]; %Matrice completa
L = Im;
n_passi = min(mA-1,nA);

%Se le righe del vettore soluzione é diverso dal numero di righe di A
if(mb ~= mA)
    error('Fattorizzazione non effettuabile.');
end


for k = 1 : n_passi
    % Controllo se elementi diagonali vicini alla soglia
    if ( abs(Ab(k,k)) <  soglia_verso_zero)
        error('Fattorizzazione non effettuabile.');
    else
        %mk = Vettore dei moltiplicatori. Vettore colonna con tante righe quante
        %quelle della matrice. Il primo elemento é zero, gli altri sono uguali
        % al rapporto tra l'elemento in posizione k e quello sulla diagonale
        mk=[zeros(k,1); Ab((k+1):end,k)/Ab(k,k)]; 
        %ekT = riga k-esima della matrice identità 
        ekT = Im(k,:);
        %Lk = matrice identità con il vettore dei moltiplicatori e tutto il
        %resto a zero
        Lk = Im - mk*ekT;
        Ab = Lk * Ab;
        invLk = Im + mk*ekT;
        L = L*invLk;
    end
end
R = triu(Ab(:,1:nA)); %Passo solo la parte di A, escludo b
y = Ab(:,end);
xs = solve_tri_sup_by_col(R,y);
end