## # Description Code Matlab
>** Aicha BOUJANDAR**


## Fonction main.m :
C'est la fonction utilisée pour obtenir les résultats présentés dans l'annexe du rapport.

Elle prend en argument le nom de l'image (par exemple 'lake'), la valeur du SNR et le patch qu'on veut utiliser pour l'algorithme SAPG sous la forme d'une liste patch (le patch sélectionné sera image(patch(1):patch(2),patch(3):patch(4)). 
Cette fonction crée l'image bruitée avec la valeur du SNR choisie, applique l'algorithme SAPG au patch sélectionné de l'image bruitée et utilise les valeurs des paramètres de régularisation obtenus pour appliquer l'algorithme du TGVdenoising à l'image entière. 

Elle retourne l'array de l'image bruitée, celui de l'image débruitée, ainsi que la liste des valeurs de paramètres le long des itérations, et de leurs moyennes à partir de l'itération N0 choisie.

## Fonction SAPG.m :
Elle applique l'algorithme SAPG à une image (ou un patch représentatif de l'image) y avec les paramètres theta1, theta2 : valeurs initiales des paramètres ,sig : l'écart type du bruit , c : le critère à utiliser (voir le jupyter notebook pour plus de détails sur les trois critèes), Nbiter : nombre d'itérations maximal ,N0 : seuil à partir duquel on calcule la moyenne des paramètres , T0 : warm-up itérations, t : thining de la chaîne de Markov provenant du prior.
Elle retourne la liste des valeurs prises par les paramètres et leurs moyennes le long des itérations.

## Fonction TGVprox.m :
Implémentation de l'opérateur proximal de la fonction TGV en Matlab. Elle retourne x = prox_TGV(y) avec theta1, theta2 les valeurs utilisées dans la fonction TGV, tau le paramètres du solveur, Nbiter le nombre d'itérations utilisées dans le solveur, et lam le paramètre de régularisation de l'opérateur proximal. Elle retourne également la valeur des fonction g1 et g2 définies dans l'article en x.

## Fonction  TGVdenoising.m :
Elle retourne la version débruitée de l'image y, avec theta1, theta2 les paramètres de régularisation, sig le bruit et Nbiter le nombre d'itérations utilisé dans le solveur.


