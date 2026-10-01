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
