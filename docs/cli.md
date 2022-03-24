# Command line tools

## cfdinfo

 <inputfile>  lit un fichier “complètement” (potentiellement tout format mais voir format plus bas) et donne quelques infos sur ce qu’il a trouvé

## ic3brief

 <inputfile> lit un fichier restart ic3 v2 ou v3 et liste toutes les sections sans lire les données : permet de savoir très rapidement ce qu’on trouve dedans (en particulier les variables)

cfdwrite_ic3v3 <inputfile>  écrit un restart ic3 v3 à partir d’un fichier qu’il a lu (inputfile) quel que soit le format (voir format) y compris ic3 v3 ; l’idée est d’y appliquer aussi des fonctions de filtre de données ou transformation, actuellement on peut effacer des variables avec --remove-cell-data LOC_LAM MU_LAM U_REL U_GRID (en fin de ligne de commande)

PS: les formats, actuellement
lecture v2 et v3 mais pas de transfert v2 to v3 ou inverse
les utilitaires génériques cfd* reconnaissent l’extension : .ic3 est recommandé pour la détection automatique des restarts. Ça se débrouille pour v2 ou v3. Si mauvaise (.out) ou pas extension, on peut imposer le format avec --fmt IC3