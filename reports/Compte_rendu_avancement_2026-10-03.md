ACTIN-ABP · 3 OCTOBRE 2026

# Compte rendu d’avancement

J’ai repris l’application à partir de la liste du 28 septembre. Le travail porte surtout sur une exploration plus claire, la cohérence des résidus entre analyses et la possibilité de retrouver l’origine des résultats.

## Les principaux changements

**Une navigation par question scientifique.** L’application est organisée en onze pages : résidus, interfaces actine–actine et ABP–actine, comparaison des sites, sensibilité mutationnelle, variants humains, propriétés chimiques et recherche structurale. Les commandes de données sont regroupées dans la documentation.

**Des correspondances corrigées.** Les analyses utilisent une numérotation commune sur la β-actine humaine P60709. Les 375 positions sont consultables, y compris sans contact observé. Les proportions par résidu sont calculées sur les PDB distinctes, pour éviter de compter plusieurs fois une même structure.

**Des vues reliées entre elles.** Les tableaux, empreintes de liaison, profils et vues 3D permettent d’examiner les contacts et de comparer les sites de plusieurs ABP. La sensibilité prédite par ProteoCast est distinguée de la conservation de séquence et des annotations cliniques des variants.

**Un jeu partagé actualisé.** Les tables retenues contiennent **2 325 interactions dans 159 PDB**, contre 2 152 interactions dans 151 PDB dans l’ancienne version publique enregistrée dans Git.

**Des résultats plus traçables.** Les calculs locaux identifiés peuvent être reconstruits automatiquement avec leurs sources et paramètres. Les données externes et les anciens fichiers dont l’origine reste incertaine sont signalés ; leur présence ne suffit plus à les considérer comme des mesures vérifiées.

## RSA : une référence humaine et une lecture simple

La RSA indique **à quel point un résidu est accessible au solvant** : une valeur élevée correspond à un résidu davantage exposé. Elle est maintenant calculée sur la **β-actine humaine 8DNH, chaîne B**, en comparant la chaîne seule et la même chaîne dans le fragment de quatre actines déposé, sans ABP.

Il s’agit d’une référence structurale précise, pas d’une moyenne de toutes nos actines. Le calcul est disponible pour 372 positions ; les positions absentes ou modifiées concernées restent sans valeur. Le choix de cette référence est à discuter selon les analyses. La méthode détaillée est conservée dans la documentation.

Comparaison Git : les dernières versions enregistrées avant le 29 septembre datent du 25 juillet 2026 pour le projet local et du 24 juillet pour le dépôt public. Le bilan ne peut donc pas isoler d’éventuels changements non enregistrés fin septembre. Plusieurs outils existaient déjà : ils ont aussi été réorganisés et corrigés.

ACTIN-ABP · SUITE

# Analyses nouvelles et suite du travail

## Comparer la disposition des actines

Pour savoir si deux actines s’assemblent comme dans un filament de référence, le contrôle **aligne la première actine puis regarde où se trouve la seconde**. Si sa position diffère, cela signale un assemblage à examiner. Un autre contrôle repère les contacts trop proches après placement d’une ABP sur un filament.

Ces outils aident à repérer les cas particuliers. Ils ne permettent pas encore, à eux seuls, de définir les interfaces « majoritaires/minoritaires » ou d’affirmer que deux partenaires se font compétition.

## Préparer et explorer les recherches FoldDisco

Le motif peut être vérifié dans une vue moléculaire, modifié puis exporté avec sa structure. Les résultats enregistrés sont séparés des nouvelles sélections. Des contrôles locaux ont été réalisés pour **206 motifs de 57 ABP**, sur un panel de 135 chaînes. Cela vérifie le fonctionnement sur ce panel ; ce n’est pas encore une recherche exhaustive ni une validation de la fonction des candidats.

## Ce qui reste à discuter ou à compléter

**1. Les références et l’interprétation des interfaces.** Confirmer l’usage de 8DNH pour la RSA et fixer les critères des interfaces majoritaires/minoritaires. Le fragment utilisé pour la RSA est fini ; la chaîne extraite conserve sa forme de filament.

**2. La couverture des résultats externes.** Les scores ProteoCast sont complets pour 49 des 57 entrées ABP contrôlées. Certains profils et certaines annotations FoldDisco restent manquants. Les dates ou paramètres originaux de quelques imports anciens n’ont pas pu être retrouvés.

**3. La validation des conclusions.** Prioriser quelques candidats FoldDisco et leurs témoins ; discuter les associations entre variants, maladies et sites de liaison avant interprétation. Les variants de signification incertaine ne sont pas automatiquement reclassés.

## Disponibilité et vérifications

Le code, la référence humaine et les résultats recalculés sont conservés sur les branches habituelles : **main** pour le projet complet et **master** pour la version partagée. Les 127 tests et l’exécution des onze pages passent localement ; le fonctionnement de l’application Cloud après redémarrage reste à confirmer.

[Dépôt de la version partagée](https://github.com/anaisdlss/actin-abp-app) · [Application Cloud](https://actin-abp-app-zqkrq5j2zqofzxgfevsmwa.streamlit.app/)

Références : liste de demandes du 28 septembre ; historique Git ; [structure humaine 8DNH](https://www.rcsb.org/structure/8DNH) (Arora et al., eLife, 2023). Les paramètres, sources et limites sont détaillés dans GUIDE.md, SCIENTIFIC_ANALYSES.md et les reçus de calcul.
