function [AAE,FFE] = elimine(AA, FF, Refneu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimine :
% Routine qui fait la pseudo-elimination pour les noeuds avec Refneu(i) = 1
% Conditions aux limites Dirichlet homogenes.
%
% SYNOPSIS elimine(AA,FF,Refneu)
%          
% INPUT * AA, FF, Refneu: La matrice et le second membre associes au
% probleme sans elimination et references des noeuds sur la frontiere
% Dirichlet.
%
% OUTPUT - AAE, FFE: On rend la matrice et le second membre une fois la 
% pseudo elimination realisee.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On copie la matrice AA et le second membre dans AAE et FFE.
AAE = AA;
FFE = FF;

% Nombre de sommets dans la discretisation.
Nbpt = size(AAE,1);

% A COMPLETER
for l=1:Nbpt
    if Refneu(l)==1  || Refneu(l) ==2  
        FFE(l)=0;
        for i=1:Nbpt
            AAE(i,l)=0;
            AAE(l,i)=0;
            if i==l
                AAE(i,l)=1;
            end
        end
    end

end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024


