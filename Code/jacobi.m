function [X] = jacobi(A, B, precision)

% A : la matrice de coefs tel que A*X = B
% B : la matrice colonne tel que A*X = B
% precision : la précision voulue sur le résultat

% Converge si A est inversible et si sa diagonale strictement dominante
iter_max = 30;      % on fixe le nb d'iterations maximale pour trouver la convergence

n = size(A, 1);     % le nombre d'equations (la taille de A)
X=rand(n, 1);       % on prend un X arbitraire

% Boucle qui recalcule X tant de fois que l'on est pas assez précis
for iterations=1:iter_max
    % On stocke le Xn précedent pour le calcul des coeffs Xn+1
    Xn=X;
    % On calcule chacun des coefficients du vecteur Xn+1
    for i=1:n
        % ATTENTION la somme se fait pour tout i != j (de 1:i-1 à i+1:n)
        X(i)=(B(i)-A(i, [1:i-1, i+1:n])*Xn([1:i-1, i+1:n]))/A(i, i);
    end    
    % si la norme 2 du vecteur (A*X-B) est inférieure à la precision alors on est assez
    % précis et on peut arrêter la boucle
    if(abs(norm(A*X-B, 2))<precision)
        iterations
        break
    end
end
