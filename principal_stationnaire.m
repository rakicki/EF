% =====================================================
% principal_stationnaire;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour :
%
% 1) l'equation suivante stationnaire, avec conditions de
% Dirichlet homogene
% | u - div(\sigma \grad u)= f,   dans \Omega=\Omega_1 U \Omega_2
% |         u = 0,   sur le bord
%
% avec
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% =====================================================

% Donnees du probleme.
nom_maillage = 'geomRect_h_0.17.msh';
affichage = true; % false; %

% Lecture du maillage et affichage.
[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri] = lecture_msh(nom_maillage);

% Declarations des matrices EF et vecteur second membre.
KK = sparse(Nbpt,Nbpt);
MM = sparse(Nbpt,Nbpt);
FF = zeros(Nbpt,1);

% Boucle d'assemblage.
for l=1:Nbtri
    % Coordonnees des sommets du triangles
  % A COMPLETER
  S1=Coorneu(Numtri(l,1),:);
  S2=Coorneu(Numtri(l,2),:);
  S3=Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l
   Kel=matK_elem(S1, S2, S3,Reftri(l));
   Mel=matM_elem(S1, S2, S3);
    for i=1:3
      I=Numtri(l,i);
      for j=1:3
        J=Numtri(l,j);
        MM(I,J)=MM(I,J)+Mel(i,j);
        KK(I,J)=KK(I,J)+Kel(i,j);
      end
    end

    
end

% Calcul du second membre.

FF=zeros(Nbpt,1);
for i=1:Nbpt
  FF(i,1)=f(Coorneu(i,1),Coorneu(i,2));
end
LL=MM*FF;
% Matrice EF complete.
AA = MM + KK;

% Pseudo-elimination.
[tilde_AA, tilde_LL] = elimine(AA,LL,Refneu);

% R?solution du prolb?me par inversion.
UU = tilde_AA \ tilde_LL;

% Visualisation de sigma et de la solution.
if affichage
    afficheSigma(Numtri, Reftri, Coorneu);
    affiche(UU, Numtri, Coorneu, sprintf('Stationnaire - %s', nom_maillage));
end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024



