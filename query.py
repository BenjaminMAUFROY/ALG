#!/usr/bin/env python3

# import
import pickle
import time
import model
import argparse
import os
import pathlib

# Arguments
parser  = argparse.ArgumentParser(
    description = "Ce programme requête un index de génomes à partir d'une ou plusieurs séquences en entrée"
)
parser.add_argument("-i", "--index", type = str, help = "Un fichier d'indexation sérialisé")
parser.add_argument("-q", "--query", type = str, help = "Un fichier fasta contenant une ou plusieurs séquences à requêter sur l'alphabet {A,C,G,T}")
parser.add_argument("-o", "--outer_file", type = str, help = "Un chemin vers un fichier en sortie de type .txt", default = "/query_result.txt")
args = parser.parse_args()

# Main
# Création d'un fichier de sortie
if args.outer_file == "/query_result.txt":
    default_path = pathlib.Path(args.query).parent
    storage_path = f"{default_path}{args.outer_file}"
else:
    if os.path.exists(args.outer_file):
        storage_path = args.outer_file
    else:
        os.mkdir(args.outer_file)
        storage_path = args.outer_file

# Désérialisation
index = pickle.load(args.index)
print(f"OUT TIME_DESERIALISATION: {time.process_time()}")

# Requêtage

