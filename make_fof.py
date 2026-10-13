#!/usr/bin/env python3

"""

make_fof.py : crée un FOF (File Of Filenames) à partir d'un chemin.

Usage :
    python make_fof.py <chemin> [-o sortie.fof] [-p MOTIF] [--no-recursive]

Exemples :
    python make_fof.py chr21_fragment.fa             -> chr21_fragment.fa.fof
    python make_fof.py ./genomes                     -> genomes.fof (tous les fichiers du dossier)
    python make_fof.py ./genomes -o fasta.fof -p "*.fa"

"""

import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Crée un FOF (file of filenames) à partir d'un chemin.")
    parser.add_argument("chemin", help="Chemin d'un fichier ou d'un dossier")
    parser.add_argument("-o", "--output", help="Fichier FOF de sortie (par défaut : <nom>.fof)")
    parser.add_argument("-p", "--pattern", default="*", help='Motif de filtrage pour un dossier, ex. "*.fa" (défaut : *)')
    parser.add_argument("--no-recursive", action="store_true", help="Ne pas descendre dans les sous-dossiers")
    args = parser.parse_args()

    chemin = Path(args.chemin).expanduser().resolve()
    if not chemin.exists():
        sys.exit(f"Erreur : '{args.chemin}' n'existe pas.")

    sortie = Path(args.output) if args.output else Path(f"{chemin.name}.fof")

    if chemin.is_dir():
        globber = chemin.glob if args.no_recursive else chemin.rglob
        fichiers = sorted(str(p) for p in globber(args.pattern) if p.is_file())
    else:
        fichiers = [str(chemin)]

    sortie.write_text("".join(f + "\n" for f in fichiers))
    print(f"FOF créé : {sortie} ({len(fichiers)} fichier(s))")


if __name__ == "__main__":
    main()