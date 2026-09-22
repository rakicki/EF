function Kel = matK_elem(S1, S2, S3, Reftri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matK_elem :
% calcul la matrices de raideur elementaire en P1 lagrange.
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%       * Reftri : reference du triangle.
%
% OUTPUT * Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul utilise une formule de quadrature de Gauss-Legendre.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Reftri ~= 1 && Reftri ~= 2
    error('Un des triangles a une réference differente de 1 et 2.');
end

% preliminaires, pour faciliter la lecture.
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps)
  error('l aire d un triangle est nulle!!!');
end
% Points et poids de quadrature.
S_hat(:,1) = [1/6; 1/6];
S_hat(:,2) = [2/3; 1/6];
S_hat(:,3) = [1/6, 2/3];
w0 = 1/6;

% Gradients des fonctions de base sur le triangle de reference.
% A COMPLETER
norm = zeros(3, 2);
norm(1, :) = [-1, -1];
norm(2, :) = [1, 0];
norm(3, :) = [0, 1];
% Transformation géométrique associée au triangle courant.
B_l = [x2-x1, x3-x1; y2-y1, y3-y1;];% A COMPLETER
S_l = [x1;y2];% A COMPLETER
Kel = zeros(3,3);
%Kel=w0*D*dot((B_l')\(norm(0,:)'),(B_l')\(norm(0,:)'));
% Boucle sur les fonctions de bases locales.
%*sigma_1(S_hat(1,q),S_hat(2,q))
for q=1:3
    Kelq=Kel;
    for i=1:3
        for j=1:3
            F_Sq=B_l*S_hat(:,q)+S_l;
            if Reftri==1     
                Kelq(i,j)=w0*sigma_1(F_Sq(1),F_Sq(2))*abs(D)*dot((B_l')\(norm(i,:)'),(B_l')\(norm(j,:)'));
            end
            if Reftri==2 
                Kelq(i,j)=w0*sigma_2(F_Sq(1),F_Sq(2))*abs(D)*dot((B_l')\norm(i,:)',(B_l')\norm(j,:)');
            end            
        end
    end
    Kel=Kel+Kelq;
end

% A COMPLETER

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024


