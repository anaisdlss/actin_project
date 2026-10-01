# Analyses complémentaires demandées par Rémi — 1 octobre 2026

Ce lot intègre les variants humains, des analyses de conservation, la régénération des figures S1 et des éléments quantitatifs et structuraux pour examiner les interfaces majoritaires/minoritaires. Les méthodes et les limites ci-dessous priment sur les notes historiques du diagnostic initial.

## Variants des six gènes

Les sources de l'application séparée sont conservées dans `data/human_variants/sources`, avec les six séquences UniProt de référence. Les dates de téléchargement initiales sont inconnues : la date d'importation ne doit pas être présentée comme une date de mise à jour de ClinVar ou gnomAD. L'import comprend 3 039 enregistrements sources, et non 3 039 variants distincts.

Les six séquences sont alignées globalement sur ACTB/P60709 (BLOSUM62 ; ouverture −10, extension −0,5). Une position n'est conservée que si sa correspondance est identique dans tous les alignements optimaux. Les résidus propres aux insertions et les discordances de lettre de référence restent dans l'audit et sont exclus des figures communes. Une observation sur une autre isoforme n'est pas assimilée à un variant clinique de P60709.

Un audit des 1 623 identifiants ClinVar compare le gène et le changement protéique à l'intitulé officiel ESummary. Il détecte quatre entrées récupérées sous ACTA2 alors que le changement concerne FAS. Ces quatre lignes sont exclues ; elles avaient également échoué au contrôle de séquence. Les anciens résultats de recherche par gène ne suffisaient donc pas à vérifier l'identité. Bilan de correspondance : 3 022 enregistrements mappés, 17 exclus pour identité/séquence/position. Le contrôle des effectifs gnomAD exclut aussi 158 lignes avec AC = 0 des analyses de présence en population. Au total, 2 864 enregistrements sont éligibles aux analyses et 175 restent uniquement dans les sources et l’audit. Les effectifs et le type exome/génome d’origine sont conservés. L'audit conserve sa date et les intitulés dans `clinvar_identity.json` ; les classifications et conditions originales restent celles de l'instantané, avec la classification actuelle conservée séparément pour comparaison.

L'application présente une heatmap par gène et catégorie, des pistes inter-gènes et une piste globale, le profil ProteoCast et le nombre de noms d'ABP contactant chaque position. Les tableaux conservent les identifiants ClinVar/gnomAD, leurs annotations et les numérotations. Les scores ProteoCast d'une substitution ne sont joints que si la lettre de référence de l'isoforme correspond à P60709 ; ils restent des scores du modèle P60709, pas des prédictions propres aux six gènes.

Les conflits explicitement annotés par ClinVar sont présentés séparément des observations inter-gènes : une substitution P/LP dans un gène et observée dans gnomAD sans entrée ClinVar dans un autre n'est pas déclarée bénigne ou tolérée. Les fractions par empreinte comptent les positions distinctes, avec des catégories pouvant se chevaucher lorsqu'elles concernent différentes substitutions. L'ASA moyenne des positions P/LP est d'abord moyennée entre chaînes d'interface, puis entre positions, uniquement parmi les contacts positifs observés.

L'association condition–empreinte utilise les mentions de conditions parmi les variants à classification agrégée P/LP d'un même gène. Le tableau 2×2 compare positions mentionnant la condition / autres positions P/LP, dans / hors empreinte. Fisher bilatéral et correction BH parmi les ABP affichées ; ces mentions agrégées ne constituent pas des assertions cliniques propres à chaque condition. Les VUS restent des VUS. [Définitions ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/) et [niveaux de revue](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).

## Conservation

Le calcul repart des 7 500 scores ProteoCast de l'actine et de sa séquence requête, retrouvés dans le projet de stage. La séquence, les lettres de référence, les 20 alternatives par position et la finitude des valeurs sont vérifiées. Le score de sensibilité est l'opposé de la moyenne des 20 scores fournis, incluant l'acide aminé inchangé à score nul, conformément à l'agrégat historique. Il ne représente ni un pourcentage d'identité ni un diagnostic clinique.

La heatmap mutationnelle, l'alignement téléchargeable, les pistes sensibilité/RSA/nombre d'ABP/contact homo, les distributions par classe de contact et les résumés par ABP et par site sont disponibles. Les corrélations de Spearman excluent les valeurs manquantes et indiquent n, p et q BH pour les quatre tests. Les résidus ne sont pas des observations indépendantes ; ces statistiques sont exploratoires et les empreintes sont influencées par l'échantillonnage structural. Les noms d'ABP sont ceux des sources, pas des familles normalisées.

Les empreintes et indicateurs de contact sont recalculés à partir des couples interaction–chaîne exacts, en incluant les deux côtés des contacts homo. M1 conserve son score ProteoCast à la colonne MAFFT 3 ; il n'est plus perdu dans une insertion. La RSA conserve son mapping structural et ses valeurs manquantes. L'origine monomère/filament de cette ancienne mesure n'est pas établie : elle ne doit pas être décrite comme une accessibilité dans le filament sans cette vérification.

## Interfaces fréquentes et contextes cofiline/coronine/WDR1

Le dénominateur de ce dépôt est 159 PDB avec contacts homo. Chaque PDB compte une seule fois par site ou cluster C70 ; les deux côtés d'une interface contribuent à leur propre site. Le contexte est défini par les noms de partenaires non-actine présents dans la PDB (cofilin, coronin, WDR1, WD-repeat-containing protein 1, actin-interacting protein 1).

| Site | PDB distinctes | Fréquence dans le jeu | Présence dans le contexte |
| --- | --- | --- | --- |
| 6685_2 | 154/159 | 96.9 % | 24/29 |
| 6685_1 | 143/159 | 89.9 % | 13/29 |
| 6685_3 | 136/159 | 85.5 % | 6/29 |
| 6685_4 | 136/159 | 85.5 % | 6/29 |
| 6685_23 | 25/159 | 15.7 % | 25/29 |
| 6685_274 | 17/159 | 10.7 % | 15/29 |
| 6685_109 | 15/159 | 9.4 % | 4/29 |
| 6685_17 | 13/159 | 8.2 % | 13/29 |

Les sites 6685_1–4 constituent bien le groupe le plus fréquent dans ces données. Le site 6685_23 et plusieurs autres sites sont associés au contexte cofiline/coronine/WDR1. Les tests et leurs listes de PDB sont exportés ; une entrée PDB n'est pas nécessairement une expérience indépendante. Une faible fréquence ne définit pas à elle seule une interface minoritaire biologique.

Un contrôle structural porte sur [3J8A, F-actine–tropomyosine](https://www.rcsb.org/structure/3J8A), [5YU8, filament décoré par la cofiline](https://www.rcsb.org/structure/5YU8) et [6VAO, filament décoré par la cofiline-1](https://www.rcsb.org/structure/6VAO). Après alignement des Cα d'une actine sur une actine de référence, le RMSD du voisin est mesuré sous la même transformation, sans le réaligner. Les deux affectations possibles des chaînes de la paire de référence sont essayées. Au moins 300 Cα communs par sous-unité sont requis.

Pour 6VAO, le RMSD médian du voisin vaut environ 0,94 Å (C70 0_7797_10) et 1,08 Å (0_7797_12) par rapport à 5YU8, contre environ 7,01 Å et 7,11 Å face aux paires correspondantes de 3J8A. Ce contrôle représentatif appuie une géométrie distincte. Il ne valide pas tous les clusters, l'état de chaque assemblage, ni des conflits stériques avec les ABP. 3J8A contient de la tropomyosine : ce n'est pas un filament nu universel.

Un bouton charge les sites observés dans 5YU8 (6685_1, 6685_2, 6685_17, 6685_23) dans la comparaison des empreintes, face à 6685_1–4. Les empreintes affichées restent l'union de toutes les observations appartenant aux sites sélectionnés ; elles ne représentent pas exclusivement cette structure. L'étiquette biologique finale « majoritaire/minoritaire » pour l'ensemble des interfaces reste à examiner avec Rémi à partir de ces résultats.

## Figures et reproductibilité

329 PNG S1 ont été régénérés dans ce dépôt : heatmaps globales, profils de chaque site, décompositions C70, empreintes homo/hétéro, sites de référence et nombre de sites ABP. L'axe couvre P60709 1–375 ; les tableaux internes conservent les coordonnées MAFFT pour ne pas casser les jointures. Les profils mono-ligne sont recalculés depuis les données corrigées, et non relus dans l'ancien CSV. Les C70 à une seule interaction sont conservés. Les deux types d'un site mixte restent séparés dans la vue globale ; le profil individuel réunit les C70 de ce site S1.

Les CSV numériques, SHA-256 des sources et du code, liste exacte des sorties et date de génération sont enregistrés dans `reports/s1_figures_all.json`. Les anciens fichiers remplacés sont conservés localement dans `reports/figure_backups` (non publiés). Les autres figures historiques C70, réseaux et exports spécialisés ne sont pas toutes régénérées par ce lot : seuls les fichiers du manifeste sont couverts.

Commandes depuis la racine du dépôt, dans l'environnement Python du projet :

```bash
python tools/audit_clinvar_identity.py
python tools/export_scientific_audit.py
python tools/check_interface_geometry.py
python tools/regenerate_scientific_figures.py
python -m unittest discover -s tests
python tools/check_app.py
```

L'audit ClinVar utilise le réseau et reprend les identifiants non encore vérifiés ; il ne rafraîchit pas automatiquement ceux déjà en cache. Les autres commandes n'actualisent pas les bases distantes. Le code est partagé entre les deux dépôts ; les figures et statistiques structurales sont recalculées sur le jeu propre à chacun. Aucun remplacement du jeu structural public par celui du projet complet n'est effectué.
