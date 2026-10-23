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
# Ouverture du fichier de requêtage
seq_list=[]
seq = ""
with open(args.query, "r") as filin:
    for line in filin:
        if line.startswith(">"):
            seq_list.append(seq)
            seq = ""
        else:
            seq += line.strip()
seq_list = seq_list[1:]

# Requêtage de l'index

kmer = 31
results = {} # Création d'un dictionnaire associant les séquences requêtes aux génomes où elles sont présentes

for i, query in enumerate(seq_list):
    results[f"query{i+1}"] = []
    for j, genome in enumerate(index):
        results[f"query{i+1}"].append("0\t")
        if genome["Filtre_Bloom"].query(query[:kmer]):
            if genome["FM_index"].backward_search(query):
                results[f"query{i+1}"][j] = "1\t"
print(f"OUT TIME_QUERY: {time.process_time()}")

# Création du fichier de sortie
i_gen = 3
with open("test.txt", "w") as querin:
    querin.write("i_gen\t")
    for i in range(i_gen):
        querin.write(f"{i+1}\t")
    querin.write("\n\n")
    for j, values in enumerate(results):
        querin.write(f"query{j+1}\t")
        for result in results[values]:
            querin.write(result)
        querin.write("\n")



    
