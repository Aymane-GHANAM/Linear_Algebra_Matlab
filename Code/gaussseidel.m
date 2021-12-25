% Fonction qui approxime la solution d'un système de n équations à n inconnues de la forme AX = B (avec A et B connues)
function [X] = gaussseidel(A, B, precision)

% A : la matrice de coefs tel que A*X = B
% B : la matrice colonne tel que A*X = B
% precision : la précision voulue sur le résultat

% Converge si A est inversible et si sa diagonale est strictement dominante
iter_max = 30;      % on fixe le nb d'iterations maximale pour trouver la convergence

n = size(A, 1);     % le nombre d'équations (la taille de A) 
X = rand(n,1);      % on prend un X arbitraire

% Boucle qui recalcule X tant que l'on a pas itéré iter_max fois
for iterations=1:iter_max
    % On calcule chacun des coefficients du vecteur Xn+1
    for i=1:n
        % Les valeurs déja calculées sont utilisées
        % Sinon on utilise les valeurs du précédent vecteur (encore dans X)
        X(i)=(B(i)-A(i, [1:i-1, i+1:n])*X([1:i-1, i+1:n]))/A(i, i);
    end
    if(abs(norm(A*X-B, 2))<precision)
        iterations
        break
    end
end

%X(i)=(B(i)-A(i, 1:i-1)*X(1:i-1)-A(i, i+1:n)*X(i+1:n))/A(i, i);
    