# Diagnostic et première réorganisation de l’application Actin–ABP

Diagnostic du 29 septembre 2026, fondé sur les deux dépôts du Bureau et sur le document de Rémi du 28 septembre 2026. Le présent lot commence la réorganisation demandée ; il ne réalise pas l’ensemble du programme scientifique du document.

## Organisation des deux dépôts

Le projet `actin_project` est la source du code et du pipeline de calcul. Le dépôt `actin-abp-app` est la version publique, avec des données préembarquées et le marqueur `data/.slim_deploy` qui désactive le pipeline. Les deux copies locales étaient propres et alignées sur leur branche distante au début du diagnostic. Leur code était identique, sauf deux corrections d’import présentes seulement dans le projet complet (`network_viz.py` et `bfactor_c70_interface.py`).

Le projet complet vient d’être cloné et ses données ne sont pas encore téléchargées, comme confirmé par Anaïs. Ce n’est pas un dysfonctionnement. La version publique contient environ 930 Mo de données, dont 2 152 lignes dans la table d’interactions et 161 sites S1 (15 homo, 144 hétéro, 2 mixtes). La table détaillée couvre 151 PDB ; le sélecteur des structures en propose 155. Cette différence entre tables doit être expliquée avant d’utiliser les totaux dans le manuscrit.

La stratégie recommandée est de développer dans le projet complet, de tester avec un jeu de données figé, puis de reporter le même code dans la version publique. Il ne faut pas développer deux implémentations scientifiques indépendantes. Une mise à jour de code est distincte d’une mise à jour de données. Ne pas remplacer les données publiques depuis un projet local vide.

## Premier lot

- Navigation en onze sections dans l’ordre demandé par Rémi, avec les analyses existantes déplacées vers leur thème : résidus, interfaces, comparaisons, conservation, propriétés et homologues.
- Gestion des téléchargements locaux sous Documentation ; le bouton de relecture du cache explique qu’il ne télécharge pas de données.
- Filtered data renommé Summary tables, avec une description pour chacun des onze tableaux et des définitions pour les principales colonnes de correspondance et d’ASA.
- Vue des sites homo et mixtes pour les interfaces actine–actine, vue des sites hétéro et mixtes pour les interfaces ABP–actine, accès au détail d’un cluster et graphe du nombre de PDB contenant des contacts homo. Aucune classification majoritaire/minoritaire n’est inférée automatiquement.
- Sélecteurs indépendants pour la conservation des ABP, la conservation de l’actine, la chimie des interfaces et FoldDisco ; tri alphabétique des principaux sélecteurs de noms et de clusters.
- Correction du tableau de comparaison par paire : les autres ABP sont recherchées dans les contacts globaux, et les deux ABP sélectionnées sont exclues de la colonne des autres partenaires.
- Message explicite lorsqu’aucun contact ABP n’est enregistré pour un résidu. Cela ne signifie pas qu’il n’existe aucun contact actine–actine.
- Protection de la construction du déploiement : vérification des fichiers sources indispensables avant toute suppression du précédent dossier deploy.
- Section Human actin variants explicitement indiquée comme non intégrée. Les analyses de variants ne sont pas présentes dans ces deux dépôts ; un autre dossier `master/M2/stage/actin_variants` existe sur le Bureau mais n’a pas été intégré ni audité dans ce lot.

La réorganisation utilise des conteneurs pour conserver l’ordre des calculs et leurs dépendances. Elle améliore l’organisation de l’écran, mais ne constitue pas encore une optimisation du temps de chargement : l’application exécute toujours ses différentes vues lors d’un rafraîchissement. Une séparation ultérieure en pages chargées à la demande est souhaitable après stabilisation des dépendances.

Les textes de GUIDE.md sont conservés, conformément à la demande de Rémi de les fournir ultérieurement. Les intitulés historiques qu’ils contiennent devront être actualisés.

## Priorité scientifique avant de nouveaux calculs

### Numérotation commune

`mafft_pipeline.build_seq_to_col_mapping` transforme une position de séquence en colonne d’alignement, pas en position de la séquence de référence sans gaps. Ajouter P60709 à un alignement n’assure donc pas à lui seul une numérotation UniProt. La propagation actuelle réutilise ces colonnes comme positions canoniques. Un repli sur une identité position→position existe aussi lorsque l’alignement manque.

Dans le jeu public examiné, l’alignement d’actine `cluster_6685.aln` ne contient pas l’entrée `P60709_ref`. D’autres alignements contiennent cette référence avec des gaps. Il faut donc auditer à la fois l’algorithme et la provenance des fichiers déjà calculés.

La table `conservation_vs_asa_per_position.csv` contient 375 lignes, avec des coordonnées allant de 5 à 380. Certaines vues limitent pourtant leur axe à 375. Supprimer simplement cette limite, ou remplacer les lettres d’acides aminés par la séquence P60709 sans remappage, risquerait de présenter des correspondances erronées.

Prochaine étape : produire et tester une table explicite PDB / chaîne / position structurale / position séquence / position P60709, avec traitement des insertions et des positions non alignables ; puis recalculer les résultats dérivés et vérifier les mêmes coordonnées dans les variants, ProteoCast et les structures 3D. Conserver un instantané de l’ancien jeu pour mesurer les différences.

### Interprétation et provenance

- RSA et ASA enfouie ne sont pas interchangeables. La méthode de calcul, la structure de référence et les unités doivent être décrites. Une valeur RSA est absente dans la table de conservation examinée.
- La fiche de résidu utilise actuellement des contacts hétéro ; la lettre affichée provient des résidus d’interface. Cela explique pourquoi une position non contactée peut ne pas avoir de lettre, sans prouver que le résidu manque dans la séquence.
- Les fréquences par PDB doivent préciser le dénominateur, dédupliquer les copies d’une même structure et prendre en compte les chaînes pertinentes. Les pourcentages signalés par Rémi restent à auditer.
- Un recouvrement de résidus ou un score Jaccard n’établit pas à lui seul un clash stérique ni une compétition fonctionnelle. Les deux niveaux d’analyse doivent rester distincts.
- Définir et faire valider les interfaces majoritaires/minoritaires et les structures de référence avant de produire leurs empreintes et des conclusions de compatibilité.
- Le pipeline saute les étapes déjà présentes : l’existence de fichiers ne garantit pas une actualisation des sources. Prévoir un manifeste daté des téléchargements, des versions d’outils, des paramètres et des résultats, ainsi qu’un mode de mise à jour explicite.

## Suite du programme de Rémi

| Lot | Travaux restants | Dépendance principale |
| --- | --- | --- |
| Données et reproductibilité | Comparer les seuils de 3 et 5 actines, qualifier les structures mutées/pathologiques et les toxines, vérifier les noms des protéines, dater les jeux de données | Critères documentés et téléchargement complet |
| Référentiel des résidus | Mapping commun, extrémités, lettres au survol, détails pour tout résidu, contacts homo et hétéro, ASA des partenaires, pourcentages par PDB | Audit de la numérotation |
| Interfaces | Seuil ASA partagé, concordance des types de clusters, empreintes majoritaires/minoritaires, RSA et recouvrements ABP/actine, distinction visuelle des partenaires | Référentiel et classification validés |
| Comparaisons | Matrice Jaccard étendue, comparaison d’empreintes et validation 3D des incompatibilités | Empreintes cohérentes |
| Conservation de l’actine | Vue générale ProteoCast, profils et corrélations avec RSA/interfaces/nombre d’ABP, détails et quantifications par cluster | Mapping ProteoCast validé |
| Variants humains | Intégration de l’application dédiée, comparaisons entre gènes, annotations conflictuelles, associations maladie/empreinte | Sources et versions des variants, mapping entre isoformes |
| Variants de signification incertaine | Analyses exploratoires et protocole de validation ; ne pas transformer automatiquement les scores en classifications cliniques | Cohortes, labels et validation indépendante |
| Conservation des ABP | Téléchargement MSA, profils complets, domaine ABD et contacts, calculs manquants et échecs | Campagne ProteoCast séparée et documentée |
| Propriétés physicochimiques | Surfaces et profils actine/partenaires, pondération des contacts par ASA à définir, comparaisons entre familles partageant un site | Définition statistique de la pondération |
| Homologues | Couverture FoldDisco, motif éditable/3D, scores intermédiaires et critères, GO, requêtes multi-clusters, validation sur ABD connus | Ressources de calcul et protocole de validation |
| Accessibilité et performance | Audit des couleurs de toutes les figures, chargement à la demande, exports et métadonnées de figures | Navigation stabilisée |

Les chimères, le deep learning et la coévolution sont des pistes exploratoires, pas des fonctionnalités réputées livrées. Les nouvelles analyses demandées dans le document sont des objectifs de travail ; elles ne constituent pas des résultats scientifiques établis.

## Vérifications du premier lot

Compilation des modules Python, tests de non-régression sur les sites mixtes, les autres partenaires et la préservation du déploiement, démarrage Streamlit avec le projet complet sans données et avec le jeu public. Trois tests de non-régression passent dans chaque dépôt. La comparaison par paire et la navigation vers le site mixte 6685_3 passent également, sans exception Streamlit. La vérification visuelle dans le navigateur intégré n’a pas pu être réalisée : son contrôle de sécurité n’était pas disponible. Les tests de l’interface ne valident pas l’exactitude de tous les calculs scientifiques existants.

Aucun pipeline de téléchargement, campagne ProteoCast ou FoldDisco n’a été lancé pour ce premier lot. Les données existantes n’ont pas été recalculées. Les modifications sont locales, sur la branche `codex/scientific-app-foundations` dans chaque dépôt ; la mise en ligne fera l’objet d’une étape distincte.

## Audit de la numérotation (résultat, 29 septembre 2026)

`script/residue_mapping_audit.py` compare, pour chaque chaîne d'actine du jeu public (3 330 chaînes), la numérotation de séquence de la table 3 avec P60709. Un décalage constant est retenu si au moins 85 % des lettres concordent (au moins 8 résidus).

- 3 249 chaînes sont alignées sur P60709 ; 81 restent non résolues (peu de résidus ou lettres discordantes). Aucune position n'est devinée pour elles.
- Décalages séquence de la chaîne → P60709 : −2 (1 055 chaînes), 0 (1 018), +1 (412), +4 (371), −1 (265), etc. La numérotation PDB n'est donc pas la numérotation UniProt.
- Après remappage, 97,1 % des lettres coïncident avec P60709.
- La colonne « canonique » MAFFT vaut la position P60709 + 4 ou + 5 (48 808 et 27 563 lignes). Cela explique les coordonnées jusqu'à 380 (375 + 5). Seulement 5,5 % des lettres coïncident si l'on lit cette colonne directement comme position UniProt.
- Sorties : `data/filtered/actin_numbering_audit_chains.csv` et `data/filtered/actin_residue_mapping_p60709.csv`. Les tables existantes et les figures ne sont pas modifiées.

Suite : décider si les vues utilisent P60709 comme numérotation affichée (soustraire 4 ou 5 ne suffit pas : utiliser la table par résidu), puis recalculer les résultats dérivés en conservant l'ancien jeu pour comparaison.

## Numérotation P60709 dans l'application (29 septembre 2026)

Les vues de l'actine affichent désormais la numérotation UniProt P60709 (1–375). La conversion est faite à l'affichage par `script/numbering.py`. Les tables internes et leurs jointures restent en colonnes MAFFT. La règle est lue dans la ligne P60709 (7pdz_I) de `data/alignments/cluster_6685.aln` : colonne 3 → M1 ; colonnes 6–236 → −4 ; 238–380 → −5. Les colonnes 1–2 et 4–5 (insertion N-terminale des actines α) et la colonne 237 n'ont pas de résidu P60709. `tests/test_numbering.py` vérifie la règle contre l'alignement et contre la table ProteoCast (375/375, hors extrémité N-terminale acide).

- Vues converties : vue d'ensemble et fiche du résidu, heatmap cliquable, comparaison par paire, variants, heatmaps S1 globales et par site, détail par position, heatmap actine × ABP, conservation de l'empreinte d'un ABP, réseaux de résidus et alignements côté actine.
- Les résidus 371–375 (colonnes 376–380), auparavant coupés par un filtre `≤ 375` appliqué aux colonnes, sont maintenant affichés.
- La lettre affichée est celle de P60709. Les acides aminés observés par organisme restent dans la fiche du résidu. Le jeu est dominé par l'actine α-squelettique de lapin ; les 12 positions où la lettre majoritaire diffère (E4, T5, T6, L16, HIC73, V129, M176, V201, I267, Y279, I287, A365) sont des différences d'isoforme connues.
- Les positions des partenaires (ABP) restent des colonnes de leur propre alignement. Elles ne sont pas des numéros UniProt ; une conversion équivalente par famille serait nécessaire pour le papier.
- Le résidu M1 n'apparaît pas dans la vue d'ensemble : l'alignement place la colonne portant la conservation ProteoCast de M1 (colonne 5) dans l'insertion N-terminale. Aucun contact n'y est enregistré.

## Bug corrigé : identifiants de chaîne comparés sans la casse

Les identifiants de chaîne PDB distinguent les majuscules des minuscules. Le code les comparait en minuscules (`.str.lower()`) dans 13 modules. Dans 10 PDB (6vec, 9q7k, 9q7l, 9q7m, 9q7n, 9y52, 9y9j, 9y9l, 9y9m, 9y9p ; 50 interactions), l'actine et son partenaire ne diffèrent que par la casse, par exemple `9y52_A` (actine) et `9y52_a` (cofiline). Des résidus du partenaire étaient alors comptés côté actine, avec des positions issues d'un autre alignement. Une chaîne partenaire pouvait aussi être classée « actine », et l'interaction comptée comme homo.

Les comparaisons sont maintenant exactes (toutes les tables utilisent le même format : préfixe PDB en minuscules, casse de la chaîne conservée). Dans les vues calculées à la volée, 7 positions parasites disparaissent. La concordance ligne à ligne des lettres d'actine avec P60709 passe de 92,5 % à 98,3 %. `tests/test_chain_case.py` couvre le cas 9y52.

À faire : les fichiers précalculés par le pipeline (PNG et CSV dérivés sous `data/filtered` et `data/visualisations`, exports `abp_site_domain`) ont été produits avec l'ancien code. Il faut les régénérer depuis le projet complet une fois les données téléchargées, puis les comparer à l'ancien jeu.


## Revue des vidéos et corrections du 1 octobre 2026

Les dix enregistrements du 29 septembre ont été examinés par captures successives des parcours. Ils montrent notamment un calcul MAFFT indisponible dans le build public, des couleurs rouge/vert dans les propriétés des acides aminés et un retour insuffisant après sélection d’un cluster. Les modifications postérieures aux vidéos déjà présentes dans les dossiers (numérotation, casse des chaînes, tableaux repliables) ont été conservées.

Ce lot ajoute :
- une confirmation du cluster sélectionné et un lien HTML vers son détail, en complément du défilement automatique ;
- des noms explicites pour les onze tableaux sources, avec leur nom de fichier toujours indiqué ;
- un histogramme des sites homo respectant réellement le tri décroissant des nombres de PDB ;
- des couleurs bleu/orange/violet pour la classification physicochimique et des graphiques plus grands ; des gradients bleus pour les partenaires des vues C70 et des interfaces 3D par ABP ;
- des commandes MAFFT désactivées dans le build public ou sans exécutable local, tout en conservant la lecture des alignements existants ; suppression du contrôle Force recompute qui pouvait déclencher un recalcul lors d’un autre changement de widget ;
- un compteur de noms de protéines non-actine utilisant l’annotation is_actin, au lieu d’exclure tout nom contenant le mot actin (ce qui excluait à tort des partenaires) ; une explication si les nombres de PDB des tables sources diffèrent.

Validation : dix tests passent sur le projet complet. Le test de casse des chaînes utilise maintenant la paire 9y52_A/9y52_a pour retrouver son identifiant, car l’interaction 1003 du jeu public est devenue 3112 dans le jeu local. Démarrage sans exception avec chacun des jeux de données ; comparaison par paire et sélection du site 6685_3 vérifiées avec le jeu public. La vérification du défilement réel dans le navigateur reste limitée par l’accès à cet outil.

Le projet complet possède désormais des données : son sélecteur présente 163 PDB et 57 noms d’ABP, contre 155 PDB et 54 noms dans le jeu public. Le code est synchronisé ; les jeux de données ne sont pas remplacés par ce lot. Il reste à vérifier la complétude des résultats dérivés avant de mettre à jour les données publiques. Pas de publication GitHub ou Cloud effectuée.

## Correspondance des résidus (1 octobre 2026)

Summary tables propose une recherche par PDB et chaîne d’actine avec export CSV : numéro structural (insertions conservées), position de séquence, acide aminé observé, colonne MAFFT et position P60709 utilisée dans l’application. Seuls les résidus d’interface présents dans les sources sont couverts. Les positions sans correspondance restent vides ; la casse des chaînes est préservée. Cette vue expose la conversion existante, sans recalculer les données ni utiliser les résultats heuristiques de l’audit des offsets.

## Contrôle direct et poursuite des demandes (1 octobre 2026)

Le navigateur intégré fonctionne : erreur d’import d’un module ancien résolue par redémarrage du serveur ; réseau de résidus rétabli en distinguant les fichiers requis et l’annotation facultative des rôles. Les 375 positions sont sélectionnables. Les RSA inconnues restent non classées dans le filtre de surface. La fiche expose les contacts actine–actine dans les deux sens et actine–ABP, les ASA des deux résidus, le type source et un export CSV. Les fréquences affichent leurs effectifs et dénominateurs. Des clés de cache ignorées par Streamlit (arguments préfixés par _) sont corrigées dans les vues concernées.

Le contrôle des heatmaps identifie une contamination entre interactions : deux listes indépendantes de chaînes et d’identifiants permettaient de retenir le mauvais côté d’un contact homo. La sélection se fait désormais sur le couple exact interaction–chaîne, dans l’app et les scripts de régénération. Le bilan avant/après sur le jeu complet est consigné dans VALIDATION.md.

Ajouts : téléchargements des alignements ABP par cluster ; FoldDisco affiche toutes les catégories par défaut et explicite ses seuils ; comparaison exploratoire ≥3/≥5 actines connectées dans Summary tables lorsque les données brutes sont présentes. Les filtres des analyses existantes ne sont pas modifiés.

Les commandes reproductibles et limites de validation sont documentées dans VALIDATION.md. Les fichiers reports/dataset_integrity.json enregistrent les sources par SHA-256. Les données Cloud restent distinctes du jeu complet ; aucune fusion dans main ni mise en production n’est effectuée ici.

Comparaison des empreintes : la section actine–actine propose deux groupes de sites configurables (référence initiale 6685_1–4), une union des positions observées sur les deux côtés des contacts homo, les résidus propres/communs, un seuil local d’ASA, les scores Jaccard contre chaque ABP et la matrice ABP × ABP téléchargeable. Aucun groupe minoritaire ni clash stérique n’est inféré automatiquement. Les tests couvrent le rattachement d’une chaîne à son site et le cas de l’union vide.


## Lot complémentaire : variants, conservation et interfaces (1 octobre 2026)

Voir [SCIENTIFIC_ANALYSES.md](SCIENTIFIC_ANALYSES.md) pour les méthodes, résultats et limites actuels. Ce lot remplace le statut antérieur « variants non intégrés » et régénère les figures S1 répertoriées dans `reports/s1_figures_all.json`. Les catégories biologiques et conclusions cliniques ne sont pas déduites automatiquement.
