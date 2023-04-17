# EDOsolvingPytohn
un simple programme python pour résoudre des equations différentielle d'ordre 1 via la méthode de Runge Kutta explicite et EulerImplicite

la méthode Runge Kutta a besoin d'une matrice de coeficient A et de deux vecteurs a,b pour fonctionner, avec bien sur l'équation differentielle
et les condiditions initiale. La matrice A DOIT être triangulaire inferieure avec des coeficients nul sur la diagonale ou les résultat seront faux.

la méthode Euler Implicite à juste besoin de l'équation differentielle, sa jacobienne selon le vecteur x et des conditions initiales.

le pas de temps h indique la précision de l'approximation.

si vous souhaitez reduire les affichages de texte lors de l'execution des simulations, vous pouvez augmenter le nombre dans la condition if juste avant l'affichage du
texte (qui vaut 500) voire lui donner une valeur négative pour ne plus avoir de texte du tout.

