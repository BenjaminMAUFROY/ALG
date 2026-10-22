#!/usr/bin/env python3

# import
import pickle
import time
import model
import argparse
import os

# Main
parser  = argparse.ArgumentParser(
    description = "Ce programme indexe les séquences nucléotidiques contenues dans un ensemble de fichiers fasta entrée et les stocks dans un fichier .txt en sortie"
)
parser.add_argument("-i", "--fof_file", type = str, help = "Un chemin vers un fichier en entrée de type File of Files", required = True)
parser.add_argument("-o", "--serialized_file", type = str, help = "Un chemin vers un fichier en sortie de type .txt", default = "./index.txt")
parser.add_argument("-m", "--indexation_method", choices = ["generalized", "bf",], help = "Méthode d'indexation. Pour plus d'informations: --help_method <generalized|bf>", required = True)
parser.add_argument("--help_methods", choices = ["generalized", "bf"], help = "Description des options de méthodes d'indexation disponibles")
args = parser.parse_args()
