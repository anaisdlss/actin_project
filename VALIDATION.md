# Vérifier l’application Actin–ABP

Une validation comporte trois niveaux : cohérence du logiciel, cohérence des données, puis validité scientifique. Un démarrage réussi ne valide pas les calculs ni les choix biologiques.

## Contrôles reproductibles

Depuis la racine du dépôt, dans l’environnement Python du projet :

```sh
python -m unittest discover -s tests -v
python tools/check_dataset.py --output reports/dataset_integrity.json
python tools/check_app.py
```

Le premier contrôle couvre notamment les chaînes sensibles à la casse, les insertions, les 375 positions P60709, les valeurs RSA inconnues, les contacts homo dans les deux sens, les fréquences par PDB et l’association exacte interaction–chaîne dans les heatmaps. Le deuxième vérifie les identifiants, les chaînes et les plages ASA et enregistre les empreintes SHA-256 des sources. Le troisième démarre l’application, parcourt M1/F375, une comparaison d’ABP et un site mixte, et vérifie l’absence d’exception. Il ne lance pas de téléchargement ni de campagne scientifique.

Si Streamlit signale une erreur d’import alors que les tests passent, arrêter puis relancer le serveur depuis l’environnement du projet. Les modules déjà importés peuvent rester anciens dans un serveur ouvert pendant une réorganisation du code. L’environnement du serveur doit aussi inclure les exécutables du projet (par exemple MAFFT), pas seulement Python.

## Vérifications dans l’interface

- Summary tables : ouvrir les correspondances PDB/chaîne, vérifier une position avec insertion et télécharger le CSV.
- Residue overview : sélectionner M1 puis F375 ; vérifier qu’une donnée absente reste absente. Ouvrir les contacts détaillés et vérifier les deux catégories actine–actine/actine–ABP. Les dénominateurs des fréquences par PDB sont explicités ; des identités différentes dans un même PDB peuvent faire dépasser 100 % à la somme des lignes.
- Comparaisons : ouvrir un site homo, un site hétéro et le site mixte 6685_3 ; vérifier le changement de cluster, la heatmap, le réseau et la vue 3D. Une capture prouve un rendu, pas l’exactitude d’une géométrie ou d’un alignement.
- Conservation : vérifier une ABP avec résultat, une sans résultat et le filtre de surface. Une RSA absente ne doit jamais être classée comme exposée.
- FoldDisco : examiner les catégories intermédiaires, les valeurs brutes et les seuils ; les catégories sont des indications exploratoires, pas une preuve de fonction biologique.

## Résultats du contrôle du 1 octobre 2026

Les 19 tests de non-régression passent dans chaque dépôt. Les six étapes du test de parcours passent avec les données propres à chaque dépôt. L’accès au navigateur intégré fonctionne désormais. Le démarrage local, le tableau de correspondance, la sélection de F375, son tableau de contacts et le réseau interactif avec structure 3D ont été observés directement.

Aucun identifiant orphelin, désaccord de chaîne de contact ou pourcentage ASA hors 0–100 n’a été trouvé dans les contrôles exécutés. Les quatre PDB 7q8b, 7q8c, 7q8s et 8oid sont présents dans la liste des structures mais exclus des interactions retenues : cluster 55649 pour les trois premiers, filtre sur le titre « fragment » pour 8oid. Le rapport distingue désormais une différence de périmètre d’un véritable détail attendu mais absent.

Le projet complet contient 2 325 interactions / 159 PDB avec détails, contre 2 152 / 151 dans le jeu public. Les listes initiales en contiennent respectivement 163 et 155. La présence d’un fichier de sortie du pipeline ne prouve ni son actualité ni sa complétude : l’interface le dit explicitement.

Une correction des heatmaps S1 associe désormais la chaîne à son interaction exacte. Sur le jeu complet, 216 cellules de profils homo, réparties sur huit sites, changent (écart maximal 46,38 points d’ASA). Les profils hétéro de cette comparaison ne changent pas. Les scripts de régénération sont corrigés également ; les anciennes images précalculées ne sont pas réputées revalidées par ce contrôle.

## Avant le papier et la mise à jour Cloud

Mise à jour du sélecteur PDB (1 octobre 2026) : son contenu est désormais construit
depuis `filtered_all_data.csv`, et non la liste du préfiltre. Les quatre entrées
7Q8B, 7Q8C, 7Q8S et 8OID n'y apparaissent donc plus ; les sources sont conservées.
Le menu compte 159 PDB dans le projet complet et 151 dans le jeu public. Tous ont
actuellement un titre et un fichier 3D disponible. Une sélection devenue absente
après une mise à jour des données est réinitialisée, avec son ancienne vue 3D.
Les titres disposent de replis sur les métadonnées locales ; leur absence seule
ne provoque pas l'exclusion d'une structure retenue.

Restent à valider : la provenance structurale de la RSA (notamment filament versus monomère), les annotations de protéines/fusions, les règles majoritaire/minoritaire, la robustesse statistique et les données variants. Les figures précalculées doivent être régénérées avec les mêmes sources et paramètres que les vues contrôlées, puis comparées aux anciennes. Les versions des outils et la date des sources doivent accompagner les figures retenues.

La comparaison de seuils peut être reproduite sur le jeu complet par :

```sh
python tools/compare_actin_thresholds.py --output reports/actin_thresholds_full_snapshot.json
```

Le filtre ≥3 retient 228 PDB contre 163 pour ≥5, et 76 contre 60 noms distincts de partenaires directs, avant les filtres ultérieurs. Ces noms ne sont pas une liste de 16 nouvelles ABP validées ; plusieurs sont des isoformes, fusions ou synonymes. Le seuil des analyses principales reste inchangé.

La comparaison d’empreintes est testée avec le site 6685_17, un seuil d’ASA de 50 % et la matrice complète activée. Elle agrège les positions présentes au-dessus du seuil dans au moins une observation ; elle ne mesure pas une fréquence biologique. Une comparaison de groupes doit consigner la liste exacte des sites et le seuil utilisé.


## Lot complémentaire : variants, conservation et interfaces (1 octobre 2026)

Voir [SCIENTIFIC_ANALYSES.md](SCIENTIFIC_ANALYSES.md) pour les méthodes, résultats et limites actuels. Ce lot remplace le statut antérieur « variants non intégrés » et régénère les figures S1 répertoriées dans `reports/s1_figures_all.json`. Les catégories biologiques et conclusions cliniques ne sont pas déduites automatiquement.

Validation du lot complémentaire : **30 tests unitaires passent** ; les **8 parcours AppTest passent** sur le jeu de ce dépôt, y compris ACTG2/VUS, conflits et empreintes vides au seuil ASA de 100 %. Inspection visuelle du navigateur sur le projet complet et de la figure S1 6685_17. Audit d’identité des 1 623 entrées ClinVar achevé. Ces contrôles ne remplacent pas la revue biologique des conclusions par les auteurs.


## Diagnostic ProteoCast (1 octobre 2026)

Dans le projet complet, la campagne de 55 ABP s'est terminée sans nouveau fichier
`4.query_ProteoCast.csv`. Les 55 fichiers d'alignement `2.ali*.fasta` téléchargés
sont vides, alors que les fichiers de structure et RSA sont présents. Le manifeste
contient aussi deux entrées sans UniProt, Cofilin et Cofilin (UNC-60B), non soumises.
Les 49 résultats déjà présents dans le jeu public sont conservés.

Un contrôle en lecture des liens MSA v6 renvoyés par l'API officielle AlphaFold pour
P23528 et Q11176 a reçu HTTP 403 (AccessDenied) pour les deux alignements. Cela
oriente vers un problème de récupération des MSA ; nous n'avons pas les journaux
internes de ProteoCast permettant d'attribuer individuellement les 55 échecs à
cette cause. La [documentation ProteoCast](https://proteocast.ijm.fr/documentation/)
indique que ce mode récupère l'alignement depuis AlphaFold. L'accès aux alignements
ou l'apport de MSA valides doit être résolu avant une nouvelle campagne.

Les diagnostics montrent désormais les fichiers incomplets et les UniProt absents.
Les prochains jobs enregistrent leur identifiant et la raison d'échec localement
(`data/proteocast/abp/.job_status/`, ignoré par Git). Le bilan et le journal de la
dernière campagne restent affichés dans la session après son rechargement. Les
compteurs distinguent le total du lot, les calculs réussis, les échecs, les tâches
en cours et celles en attente. Les erreurs de préparation et les codes de sortie
ne sont plus présentés comme une réussite. Un fichier de scores vide ne compte
pas comme un résultat disponible.

Validation de cette correction : 7 tests unitaires ProteoCast réussis, 7 scénarios
AppTest isolés avec soumissions simulées réussis, et démarrage sans exception des
deux applications avec leurs données respectives. Aucun nouveau calcul externe
n'a été soumis pendant ce diagnostic. Ces vérifications valident le suivi logiciel,
pas le fonctionnement actuel du calcul sur le serveur distant.

## Validation du complément interfaces, variants et conservation — 1 octobre 2026

Les contrôles sont exécutés séparément sur les deux jeux de données. **77 tests unitaires passent dans chaque dépôt**. Les **12 parcours AppTest** de `tools/check_app.py` passent dans chaque dépôt, dont les nouveaux filtres RCSB, la comparaison inter-gènes, l'empreinte alignée Cofilin-1, les surfaces de chimie générées et les seuils sans contacts. L'audit d'intégrité ne détecte aucun échec ; seul persiste le signalement des quatre structures hors périmètre d'interactions, déjà expliqué ci-dessus.

Les smokes dédiés vérifient aussi les profils ABP (Cofilin-1, Filamin-A, filtre RSA et résultats absents), tous les onglets de contacts des sites 6685_0/6685_3/6685_2 et l'analyse Cofilin, les filtres d'annotations, les requêtes FoldDisco éditées valides/invalides/sans résultat et les seuils de RSA de référence. Le navigateur a permis d'inspecter réellement les surfaces de chimie et le profil public Cofilin-1 avec contacts et domaines ; AppTest seul ne valide pas WebGL.

L'ancien calcul de contacts MSA sélectionnait des chaînes sans restreindre leurs interactions au site choisi, et privilégiait une orientation. Après correction sur le jeu complet, le périmètre de 6685_0 passe de 17 623 à 10 815 lignes orientées ; celui de 6685_3 de 23 987 à 52 071, et de 6685_2 de 58 098 à 26 206. Ces écarts incluent la récupération des côtés homo manquants. Les aires censurées, absentes, non finies ou non positives ne sont plus inventées. Les profils pondérés utilisent des moyennes par chaîne puis PDB observé ; la composition garde les contributions de chaque acide aminé, pas uniquement l'AA dominant. Le panneau AA-pairs indique explicitement son agrégation d'aires poolées ; les vues de présence conservent aussi les contacts sans aire mesurable. Un type de contact vide n'est plus interprété comme une liaison de van der Waals.

Dans la nouvelle vue de chimie, une relecture a également imposé une moyenne par interface avant PDB puis ABP. Par rapport au premier calcul provisoire qui regroupait directement les paires d'une PDB, 234 positions changent, avec un écart maximal de fraction de 0,175746 et trois changements de classe dominante. Le code final et ses tests incluent cette correction. La projection de 7PDZ I contrôle 371 résidus et 2 894 identifiants atomiques distincts ; HIC73 est distinguée de l'histidine standard.

Les méthodes, paramètres, couvertures et limites figurent dans SCIENTIFIC_ANALYSES.md. Le nouveau calcul d'accessibilité a son propre manifeste et son contrôle numérique 480/960 points. Les tableaux de comparaison inter-gènes et de chimie sont reproductibles avec `tools/export_scientific_audit.py`. Les anciens exports spécialisés de chimie/FoldDisco non listés dans ces manifestes ne sont pas réputés régénérés. Les tests ne constituent pas une validation biologique de toutes les interfaces ou une validation clinique des variants.

## Sélections interactives et échelles de couleur — 1 octobre 2026

Les clics Plotly sont traités par des callbacks d’événement : une ancienne sélection ne réécrit plus les choix manuels des résidus ou des clusters lors des reruns. La vue d’ensemble garde le résidu choisi dans le menu, la fiche, la surface et son repère sur le graphique.

La heatmap S1 affiche une légende colorée avant le grand graphique et une échelle native compacte conservée dans les exports. Le mode absolu respecte désormais la plage annoncée 0–100 % ; le mode relatif utilise 0–1 par cluster. Les valeurs scientifiques et les pourcentages absolus au survol sont inchangés.

La matrice des variants individuels est explicitement binaire (substitution enregistrée ou absente du snapshot sélectionné). Elle n’affiche plus de dégradé continu trompeur entre 0 et 1. Le graphique agrégé utilise un vrai dégradé du nombre entier de substitutions distinctes par position ; il ne représente ni la gravité, ni une fréquence allélique.

Validation : **88 tests passent dans chacun des deux dépôts**, dont les événements Plotly transmis au moteur Streamlit, puis les choix manuels et reruns indépendants. Le smoke de l’application réelle passe au démarrage, pour le résidu P60709 47 et pour ACTA1/pathogenic. Vérification navigateur du passage clic G46 → choix manuel M47 et de la légende S1 visible.

## Navigation par pages — 1 octobre 2026

Les onze rubriques du document de Rémi deviennent onze pages natives. Une vue
sélectionnée charge ses analyses ; les vues des autres rubriques ne sont plus
construites à chaque interaction. Les calculs de données restent dans
Documentation, tandis que les analyses restent accessibles par un sélecteur
visible en haut de chaque page. Les exports et réseaux globaux sont regroupés.

Les choix explicitement identifiés de résidu, protéine, cluster et paramètres
sont conservés dans la session. Les événements Plotly et boutons de calcul ne
sont pas conservés comme des choix : revenir à une page ne rejoue pas un calcul.
Une heatmap ouvre le cluster sélectionné ; un nœud ABP ouvre la page du partenaire.

`tools/check_app.py` parcourt les onze pages et toutes leurs vues, puis vérifie
les changements de paramètres, les transitions cluster/protéine, le clic sur
une heatmap, les contrôles de calcul réservés au projet complet et le retour
au résidu M47. Il utilise les données installées ; seul le téléchargement de
secours AlphaFold est neutralisé dans ce parcours pour permettre le contrôle
hors ligne. Les résultats locaux et les autres analyses ne sont pas simulés.
Les 88 tests de régression passent dans les deux dépôts.

Les formules, seuils et jeux de données ne sont pas modifiés par cette
réorganisation. Une meilleure navigation ne remplace pas les validations
biologiques encore répertoriées dans SCIENTIFIC_ANALYSES.md.

Contrôle final : les parcours complets passent avec les jeux local et public,
y compris les réseaux de compétition et de coprésence. Dans le navigateur,
le résidu M47 est conservé après aller-retour par Summary tables ; le bouton
d’exploration du jeu public ouvre bien la page du cluster sélectionné.

## Réduire les niveaux de navigation — 1 octobre 2026

Documentation affiche directement le guide, avec un seul panneau secondaire
Data management pour le cache, le pipeline et les commandes ProteoCast. Le
sélecteur Guide / Data and calculations et le panneau de mise à jour imbriqué
sont supprimés. Les tableaux homo/hétéro, le réseau ABP sélectionné et les
résultats Structural evidence s'affichent directement. Les onze rubriques et
les choix entre analyses distinctes sont conservés. Ce lot modifie uniquement
la présentation, sans changer les données ni les formules.

## Suivi des demandes de Rémi — 2 octobre 2026

Les détails des sites sont accessibles dans leurs pages actine–actine ou
ABP–actine. La comparaison des empreintes s'affiche directement ; la
conservation d'un résidu se consulte dans la page Conservation en gardant la
sélection. La nouvelle surface homo/hétéro/mixte utilise les contacts à ASA
positive et un mapping de séquence P60709. Une seule surface est colorée par
atome ; les positions sans observation ou non mappées restent grises.

Validation finale sur les deux dépôts : **95 tests unitaires réussis** chacun,
puis **toutes les pages, leurs vues et les interactions de `tools/check_app.py`
réussies** avec les jeux respectifs. Le test ciblé `tools/check_followup.py`
vérifie l'affichage des 2 000 alignements sauvegardés pour Inverted formin-2 /
6685_213, la combinaison de sites sur 9AZ4 G (41 positions, limite serveur
signalée), et l'ouverture de la conservation de M47 sans perte de sélection.
Aucune soumission distante n'est lancée par ces tests.

Les empreintes SHA-256 de toutes les requêtes sauvegardées sont conformes aux
coordonnées exportées ; chaque position du motif possède un C-alpha exporté.
Les rapports géométriques correspondent aux tableaux sources actuels :
1 418 paires / 23 clusters dans le complet, 1 340 / 21 dans le public.
Les comptes ProteoCast et les résultats des nouvelles recherches FoldDisco
sont détaillés dans SCIENTIFIC_ANALYSES.md et les rapports d'audit datés.
Un résultat vide, un échec et un motif non pris en charge restent distingués.

Le contrôle visuel du navigateur n'a pas pu être refait pour ce lot : l'outil
a refusé l'accès faute de pouvoir vérifier une règle de sécurité administrateur.
Ce contrôle n'a pas été contourné. AppTest vérifie l'exécution et les interactions
Streamlit, pas le rendu JavaScript effectif de la nouvelle surface 3D. Une
inspection visuelle de cette surface reste donc à faire.

Reproduire les contrôles depuis chaque dépôt, avec l'environnement Python du projet :

```bash
python -m unittest discover -s tests
python tools/check_app.py
python tools/check_followup.py
```

Les nouveaux calculs ProteoCast restent bloqués par des alignements distants
vides, et deux entrées du manifeste complet n'ont pas d'accession. La validation
biologique des interfaces et des candidats ABD reste du ressort des auteurs ;
les mesures et leurs limites sont disponibles pour cette revue. Aucun
redéploiement de l'application Cloud n'est réalisé par ce lot.


## Display, terminology and 3D provenance — 2 October 2026

- 102 unit tests pass in each repository, including numeric interaction order,
  per-cell hover colours without changing measurements, real py3Dmol element
  resizing, complete-score/sequence requirements and cache refresh after import.
- `tools/check_app.py` passes across the 11 pages and their representative views
  in both datasets. `tools/check_followup.py` passes for saved FoldDisco results,
  coherent combined motifs, server-size limits and residue navigation.
- `tools/check_display.py` passes in both datasets: numeric source table, individual
  Myosin-6 pair in 6685_6, aligned partners, explicit 3D metric legend, visible zoom
  controls and cell-colour tooltip payloads. Dataset-specific ABP availability is
  preserved; the public test does not assume Actin-interacting protein 1 is present.
- The complete dataset supplies seven observed partner pairs for 6685_6. Their
  superposition uses sequence-matched actin Cα atoms; colours in a single-pair view
  use that pair's ASA measurements, rather than pooled cluster measurements.
- The existing Actin-interacting protein 1 PDB is an original AlphaFold file. Its
  B-factors are now labelled pLDDT confidence, rather than ProteoCast sensitivity.
- `git diff --check` passes. No clinical classifications or source contact tables
  were changed by this display update.

Visual browser inspection was blocked because the browser tool could not verify
its required security policy. No alternate browser path was used. AppTest validates
page execution and generated controls, but does not certify the rendered WebGL
framing or tooltip appearance on the user's screen. These visual changes therefore
still need a browser check once that tool access is available.
