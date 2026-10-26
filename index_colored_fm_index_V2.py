#!/usr/bin/env python3
"""Indexation de génomes sous forme de FM-index et sérialisation de l'index.

Ce script lit un fichier de type File of Files (FoF), c'est-à-dire un fichier
listant les fichiers FASTA à indexer, construit un index à l'aide du module
local ``model``, puis sérialise cet index sur disque avec ``pickle``.

Deux méthodes d'indexation sont prévues (option ``-m`` / ``--indexation_method``) :

- ``generalized`` : les génomes sont concaténés puis indexés ensemble sous la
  forme d'un FM-index unique accompagné d'un color_array. Cette méthode n'est
  pas encore implémentée.
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
    L'index construit, sérialisé avec ``pickle`` (le fichier est binaire malgré
    l'extension .txt), ainsi que les temps CPU de construction et de
    sérialisation affichés sur la sortie standard.

Dépendances :
    Le module local ``model``, qui doit fournir la fonction ``bf``.
"""

# Imports
import argparse     # Gestion des arguments en ligne de commande
import os           # Vérification et création de chemins
import pathlib      # Manipulation de chemins (récupération du dossier parent)
import pickle       # Sérialisation de l'index sur disque
import sys          # Arrêt du programme avec un code de sortie
import time         # Mesure du temps d'exécution

import model        # Module local contenant les méthodes d'indexation (ex. model.bf)

# Arguments
# Définition de l'interface en ligne de commande.
parser  = argparse.ArgumentParser(
    description = "Ce programme indexe les séquences nucléotidiques contenues dans un ensemble de fichiers fasta entrée et les stocks dans un fichier .txt en sortie"
)
# Fichier d'entrée : File of Files listant les fichiers fasta à indexer
parser.add_argument("-i", "--fof_file", type = str, help = "Un chemin vers un fichier en entrée de type File of Files", metavar = "filin")
# Fichier de sortie : la valeur par défaut "/index.txt" sert de marqueur pour
# créer le fichier dans le dossier du FoF (voir la partie Sérialisation)
parser.add_argument("-o", "--serialized_file", type = str, help = "Un chemin vers un fichier en sortie de type .txt", default = "/index.txt", metavar = "filout")
# Méthode d'indexation à utiliser
parser.add_argument("-m", "--indexation_method", choices = ["generalized", "bf",], help = "Méthode d'indexation. Pour plus d'informations: --help_methods <generalized|bf>")
# Option d'aide : affiche la description d'une méthode au lieu de lancer l'indexation
parser.add_argument("--help_methods", choices = ["generalized", "bf"], default = "", help = "Description des options de méthodes d'indexation disponibles")
# Lecture et analyse des arguments passés en ligne de commande
args = parser.parse_args()

# Description des méthodes
# Si --help_methods est renseigné, on affiche la description de la méthode
# choisie puis on quitte : aucune indexation ni sérialisation n'est lancée.
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
# (ils ne sont pas "required" dans argparse pour que --help_methods fonctionne seul)
if args.fof_file is None or args.indexation_method is None:
    parser.error("les options -i/--fof_file et -m/--indexation_method sont obligatoires (sauf avec --help_methods).")

# Main
# Indexation

# Indexation avec la méthode "generalized"
# (concaténation des génomes, un seul FM-index et un color_array)
if args.indexation_method == "generalized":
    # Méthode pas encore implémentée : on quitte avec un code d'erreur
    sys.exit("A faire : la méthode 'generalized' n'est pas encore implémentée.")

# Indexation avec la méthode "bf"
# (un FM-index et un filtre de Bloom par génome)
elif args.indexation_method == "bf":
    # Temps CPU du processus avant la construction
    start_build = time.process_time()
    # Construction de l'index à partir du fichier FoF
    index = model.bf(args.fof_file)
    # Durée CPU de la construction de l'index
    print(f"OUT TIME_BUILD: {time.process_time() - start_build}")

# Serialisation
# Vérification de la présence de l'option -o. Si non, création d'un fichier par défaut dans le dossier contenant le FoF.
if args.serialized_file == "/index.txt":
    # Option -o non renseignée : on utilise le dossier contenant le FoF
    default_path = pathlib.Path(args.fof_file).parent
    # Chemin final : <dossier_du_FoF>/index.txt
    storage_path = f"{default_path}{args.serialized_file}"
else:
    # Option -o renseignée : vérification de l'existence du fichier de sortie
    if os.path.exists(args.serialized_file):
        # Le fichier existe déjà : on le réutilise tel quel
        storage_path = args.serialized_file
    else:
        # Le fichier n'existe pas : création des dossiers parents manquants avec os.makedirs.
        # On passe le dossier parent et non le chemin complet, sinon un dossier
        # portant le nom du fichier serait créé à la place du fichier.
        parent_dir = os.path.dirname(args.serialized_file)
        # parent_dir est vide si seul un nom de fichier est donné (ex. : out.txt)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok = True)
        # Le fichier lui-même est créé à l'ouverture en écriture, lors de la sérialisation
        storage_path = args.serialized_file

# Création du fichier de sortie avec pickle
# Le fichier est ouvert en écriture binaire ("wb") : pickle.dump attend un
# objet fichier et non un chemin. Un fichier existant est écrasé.
start_serialisation = time.process_time()
with open(storage_path, "wb") as file_out:
    pickle.dump(index, file_out)
# Durée CPU de la sérialisation
print(f"OUT TIME_SERIALISATION: {time.process_time() - start_serialisation}")
