#!/usr/bin/env python3

# import
import pickle
import time
import model
import argparse
import os
import pathlib

# Main
parser  = argparse.ArgumentParser(
    description = "Ce programme indexe les séquences nucléotidiques contenues dans un ensemble de fichiers fasta entrée et les stocks dans un fichier .txt en sortie"
)
parser.add_argument("-i", "--fof_file", help = "Un chemin vers un fichier en entrée de type File of Files", metavar = "filin")
parser.add_argument("-o", "--serialized_file", type = str, help = "Un chemin vers un fichier en sortie de type .txt", default = "/index.txt", metavar = "filout")
parser.add_argument("-m", "--indexation_method", choices = ["generalized", "bf",], help = "Méthode d'indexation. Pour plus d'informations: --help_method <generalized|bf>")
parser.add_argument("--help_methods", choices = ["generalized", "bf"], default = "", help = "Description des options de méthodes d'indexation disponibles")
args = parser.parse_args()

# Description des méthodes
if args.help_methods == "generalized":
    print(
        "Indexation sous forme de FM_index et d'un color_array d'un ensemble de génomes par concaténation de ceux-ci."
    )
elif args.help_methods == "bf":
    print(
        "Chaque génome est indexé sous la forme d'un FM_index et d'un filtre de Bloom associé"
    )


# Indexation 
elif args.help_methods == "":
    
    # Indexation avec la méthode "generalized"
    if args.m == "generalized":
        print("A faire")

    # Indexation avec la méthode "bf"
    elif args.m == "bf":
        index = model.bf(args.filin)
        print(time.process_time())

# Serialisation
# Vérification de l'existence du fichier de sortie. Si non, création d'un fichier par défaut dans le dossier contenant le FoF.
if args.filout == "/index.txt":
    default_path = pathlib.Path(args.filin).parent
    storage_path = f"{default_path}{args.filout}"
else:
    if os.path.exists(args.filout):
        storage_path = args.filout
    else:
        os.mkdir(args.filout)
        storage_path = args.filout
        
pickle.dump(index, storage_path)
print(time.process_time())
