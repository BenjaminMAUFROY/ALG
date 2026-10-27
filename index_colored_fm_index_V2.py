#!/usr/bin/env python3
"""Indexation de génomes sous forme de FM-index et sérialisation de l'index.

Ce script lit un fichier de type File of Files (FoF), c'est-à-dire un fichier
listant des chemins vers des fichiers FASTA à indexer, construit un index à
l'aide du module local ``model``, puis sérialise cet index sur disque avec
``pickle``.

Deux méthodes d'indexation sont prévues (option ``-m`` / ``--indexation_method``) :

- ``generalized`` : les génomes sont concaténés puis indexés ensemble sous la
  forme d'un FM-index unique accompagné d'un color_array.
- ``bf`` : chaque génome est indexé séparément sous la forme d'un FM-index
  associé à un filtre de Bloom.

Utilisation :
    python3 index_colored_fm_index.py -i <fichier_fof> -m <generalized|bf> [-o <fichier_sortie>]
    python3 index_colored_fm_index.py --help_methods <generalized|bf>

Arguments :
    -i, --fof_file            Chemin vers le fichier FoF en entrée.
    -o, --serialized_file     Chemin vers le fichier de sortie. S'il n'existe
                              pas, il est créé (ainsi que les dossiers parents
                              manquants). Si l'option n'est pas renseignée,
                              un fichier ``index.txt`` est créé dans le dossier
                              contenant le FoF.
    -m, --indexation_method   Méthode d'indexation : ``generalized`` ou ``bf``.
    --help_methods            Affiche la description d'une méthode d'indexation
                              puis quitte sans effectuer d'indexation.

Sortie :
    L'index construit, sérialisé avec ``pickle``.
    Le temps CPU de construction de l'index affiché dans la console.
    Le temps de sérialisation affiché dans la console.

Dépendances :
    Le module local ``model``, qui doit fournir la fonction ``bf``.
"""

# Imports
import argparse
import model     
import os           
import pathlib
import pickle
import sys
import time

# Arguments
# Définition de l'interface en ligne de commande.
parser  = argparse.ArgumentParser(description = "Ce programme indexe les séquences nucléotidiques contenues dans un ensemble de fichiers fasta entrée et les stocks dans un fichier en sortie")
parser.add_argument("-i", "--fof_file", type = str, help = "Un chemin vers un fichier en entrée de type File of Files", metavar = "filin")
parser.add_argument("-o", "--serialized_file", type = str, help = "Un chemin vers un fichier .pkl en sortie", default = "/index.pkl", metavar = "filout")
parser.add_argument("-m", "--indexation_method", choices = ["generalized", "bf",], help = "Méthode d'indexation. Pour plus d'informations: --help_methods <generalized|bf>")
parser.add_argument("--help_methods", choices = ["generalized", "bf"], default = "", help = "Description des options de méthodes d'indexation disponibles")
args = parser.parse_args()

# Description des méthodes
# Si --help_methods est renseigné aucune indexation ni sérialisation n'est
# lancée.
if args.help_methods == "generalized":
    print(
        "Indexation sous forme de FM_index et d'un color_array d'un ensemble de génomes par concaténation de ceux-ci."
    )
    sys.exit(0)

elif args.help_methods == "bf":
    print(
        "Chaque génome est indexé sous la forme d'un FM_index et d'un filtre de Bloom associé"
    )
    sys.exit(0)

# Vérification des arguments obligatoires pour l'indexation
if args.fof_file is None or args.indexation_method is None:
    parser.error("les options -i/--fof_file et -m/--indexation_method sont obligatoires (sauf avec --help_methods).")

# Main
# Indexation

# Indexation avec la méthode "generalized"
# (concaténation des génomes, un seul FM-index et un color_array)
if args.indexation_method == "generalized":
    start_build = time.process_time()

    # Concaténation des génomes du FoF
    sequence, genome_of, paths = model.build_from_fof(args.fof_file)

    # FM-index + color array
    fm = model.FM_Index(sequence, genome_of)
    index = {"FM_index": fm.fm_index(), "Color_array": fm.get_color(), "paths": paths,}

    # Durée CPU de la construction de l'index
    print(f"OUT TIME_BUILD: {time.process_time() - start_build}")

# Indexation avec la méthode "bf"
# (un FM-index et un filtre de Bloom par génome)
elif args.indexation_method == "bf":
    start_build = time.process_time()

    # Construction de l'index à partir du fichier FoF
    index = model.bf(args.fof_file)

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
            os.makedirs(parent_dir, exist_ok = True)

        storage_path = args.serialized_file

# Création du fichier de sortie avec pickle
start_serialisation = time.process_time()

with open(storage_path, "wb") as file_out:
    pickle.dump(index, file_out)

# Durée CPU de la sérialisation
print(f"OUT TIME_SERIALISATION: {time.process_time() - start_serialisation}")
