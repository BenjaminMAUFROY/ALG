"""Indexation de génomes sous forme de FM-index et sérialisation de l'index.

Ce script lit un fichier de type File of Files (FoF), c'est-à-dire un fichier
listant des chemins vers des fichiers FASTA à indexer, construit un index à
l'aide des classes et fonctions de ``model``, puis sérialise cet index sur
disque avec ``pickle``.

Deux méthodes d'indexation sont prévues (option ``-m`` /
``--indexation_method``) :

- ``generalized`` : les génomes sont concaténés puis indexés ensemble sous la
  forme d'un FM-index unique accompagné d'un color_array.
- ``bf`` : chaque génome est indexé séparément sous la forme d'un FM-index
  associé à un filtre de Bloom.

Dans les deux cas, les séquences sont restreintes aux bases A, C, G et T.

Paramètres du filtre de Bloom (méthode ``bf``) :
    Ils se règlent en tête du script, dans trois variables :

    - ``n_hashage`` : nombre de fonctions de hachage H (7 par défaut) ;
    - ``fp_rate``   : taux de faux positifs visé (0.01 par défaut, soit 1 %);
    - ``kmer_size`` : taille des kmers insérés dans le filtre (31 par
                      défaut).

    ATTENTION : ``n_hashage`` et ``kmer_size`` doivent être identiques aux
    variables ``H`` et ``kmer`` de ``query.py``.

Utilisation :
    python3 index_colored_fm_index.py -i <fichier_fof> -m <generalized|bf>
        [-o <fichier_sortie>]
    python3 index_colored_fm_index.py --help_methods <generalized|bf>

Arguments :
    -i, --fof_file            Chemin vers le fichier FoF en entrée.
    -o, --serialized_file     Chemin vers le fichier de sortie. S'il n'existe
                              pas, il est créé (ainsi que les dossiers
                              parents manquants). Si l'option n'est pas
                              renseignée, un fichier ``index.pkl`` est créé
                              dans le dossier contenant le FoF.
    -m, --indexation_method   Méthode d'indexation : ``generalized`` ou
                              ``bf``.
    --help_methods            Affiche la description d'une méthode
                              d'indexation puis quitte sans effectuer
                              d'indexation.

Sortie :
    L'index construit, sérialisé avec ``pickle``.
    Le temps CPU de construction de l'index affiché dans la console.
    Le temps de sérialisation affiché dans la console.

Dépendances :
    Le module local ``model``, qui doit fournir la fonction ``bf`` (méthode
    ``bf``), ainsi que la fonction ``build_from_fof`` et la classe
    ``FM_Index`` (méthode ``generalized``).
"""

# Imports
import argparse
import os
import pathlib
import pickle
import time

import model

# Paramètres du filtre de Bloom (méthode "bf"), à modifier ici.
# n_hashage et kmer_size doivent être identiques à H et kmer dans query.py.

# Nombre de fonctions de hachage H
n_hashage = 7
# Taux de faux positifs visé, en fraction (0.01 = 1 %)
fp_rate = 0.01
# Taille des kmers insérés dans le filtre
kmer_size = 31

# Arguments
# Définition de l'interface en ligne de commande.
parser = argparse.ArgumentParser(
    description=(
        "Ce programme indexe les séquences nucléotidiques contenues dans un "
        "ensemble de fichiers fasta entrée et les stocks dans un fichier en "
        "sortie"
    )
)
parser.add_argument(
    "-i", "--fof_file", 
    type=str,
    help="Un chemin vers un fichier en entrée de type File of Files",
    metavar="filin",
    required=True
)
parser.add_argument(
    "-o", "--serialized_file", 
    type=str,
    help="Un chemin vers un fichier .pkl en sortie",
    default="/index.pkl", 
    metavar="filout",
    required=True
)
parser.add_argument(
    "-m", "--indexation_method",
    choices=["generalized", "bf"],
    help=("Méthode d'indexation. Pour plus d'informations: "
          "--help_methods <generalized|bf>"),
    required=True
)
parser.add_argument(
    "--help_methods", choices=["generalized", "bf"],
    default="",
    help="Description des options de méthodes d'indexation disponibles"
)
args = parser.parse_args()

# Description des méthodes
if args.help_methods == "generalized":
    print(
        "Indexation sous forme de FM_index et d'un color_array d'un ensemble "
        "de génomes par concaténation de ceux-ci."
    )

elif args.help_methods == "bf":
    print(
        "Chaque génome est indexé sous la forme d'un FM_index et d'un filtre "
        "de Bloom associé"
    )

# Indexation

# Indexation avec la méthode "generalized"
# (concaténation des génomes, un seul FM-index et un color_array)
if args.indexation_method == "generalized":
    start_build = time.process_time()

    # Concaténation des génomes du FoF
    sequence, genome_of, paths = model.build_from_fof(args.fof_file)

    # FM-index + color array
    fm = model.FM_Index(sequence, genome_of)
    index = {"FM_index": fm.fm_index(), "Color_array": fm.get_color(),
             "paths": paths}

    # Durée CPU de la construction de l'index
    print(f"OUT TIME_BUILD: {time.process_time() - start_build}")

# Indexation avec la méthode "bf"
# (un FM-index et un filtre de Bloom par génome)
elif args.indexation_method == "bf":
    start_build = time.process_time()

    # Construction de l'index à partir du fichier FoF, avec les paramètres du
    # filtre de Bloom définis en tête de script
    index = model.bf(args.fof_file, h=n_hashage, fp_rate=fp_rate,
                     kmer=kmer_size)

    # Durée CPU de la construction de l'index
    print(f"OUT TIME_BUILD: {time.process_time() - start_build}")

# Serialisation
# Vérification de la présence de l'option -o. Si non, création d'un fichier
# par défaut dans le dossier contenant le FoF.
if args.serialized_file == "/index.pkl":
    default_path = pathlib.Path(args.fof_file).parent
    storage_path = f"{default_path}{args.serialized_file}"

else:
    # Si -o renseignée : vérification et gestion de l'existence
    # du fichier de sortie
    if os.path.exists(args.serialized_file):
        storage_path = args.serialized_file

    else:
        parent_dir = os.path.dirname(args.serialized_file)

        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)

        storage_path = args.serialized_file

# Création du fichier de sortie avec pickle
start_serialisation = time.process_time()

with open(storage_path, "wb") as file_out:
    pickle.dump(index, file_out)

# Durée CPU de la sérialisation
print(f"OUT TIME_SERIALISATION: {time.process_time() - start_serialisation}")
