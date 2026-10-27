#!/usr/bin/env python3

# import
import pickle
import time
import argparse
import model
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
start_deserialisation = time.process_time()
with open(args.index, "rb") as filin:
    index = pickle.load(filin)

print(f"OUT TIME_DESERIALISATION: {time.process_time() - start_deserialisation}")

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
    if seq :
        seq_list.append(seq)
seq_list = seq_list[1:]

# Requêtage de l'index

# Détection de la méthode d'indexation : pour generalized = dictionnaire
generalized_index = None
if isinstance(index, dict) and "Color_array" in index:
    generalized_index = index
    index = []  # la boucle des filtres de Bloom ci-dessous n'a alors rien à parcourir



kmer = 31
results = {} # Création d'un dictionnaire associant les séquences requêtes aux génomes où elles sont présentes

start_query = time.process_time()

for i, query in enumerate(seq_list):
    results[f"query{i+1}"] = []
    for j, genome in enumerate(index):
        results[f"query{i+1}"].append("0\t")
        if model.query(genome["Filtre_Bloom"], query[:kmer]):
            if model.backward_search(genome["FM_index"], query):
                results[f"query{i+1}"][j] = "1\t"


# Requêtage pour index generalized
if generalized_index is not None:
    fm_index = generalized_index["FM_index"]
    color = generalized_index["Color_array"]
    nb_genomes = len(generalized_index["paths"])

    for i, query in enumerate(seq_list):
        results[f"query{i+1}"] = ["0\t"] * nb_genomes
        interval = model.backward_search_interval(fm_index, query)

        if interval is not None:
            start, stop = interval
            for genome_index in set(int(c) for c in color[start:stop + 1] if c is not None):
                results[f"query{i+1}"][genome_index] = "1\t"


print(f"OUT TIME_QUERY: {time.process_time() - start_query}")
print(results)

# Création du fichier de sortie

i_gen = len(index)

with open(args.outer_file, "w") as querin:
    querin.write("i_gen\t\t")
    for i in range(i_gen):
        querin.write(f"{i+1}\t")
    querin.write("\n\n")
    for j, values in enumerate(results):
        querin.write(f"query{j+1}\t\t")
        for result in results[values]:
            querin.write(result)
        querin.write("\n")



    
