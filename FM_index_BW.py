"""
FM-index généralisé via Burrows-Wheeler (BWT) et color-array.

Approche 1:
N génomes concaténés --> BWT --> FM index + Color array

On étend la notion de FM-index aux collections de séquences en construisant le FM-index (BW T, Count, Rank)1
de la chaı̂ne de caractères S obtenue en concaténant les N génomes, séparés par le caractère spécial $
(également présent en fin de chaı̂ne). On associe à ce FM-index, un color array. C’est une table, similaire
au suffix array, mais au lieu de stocker la position de chaque suffixe, on stocke l’identifiant (la couleur) du
génome (un entier dans [0, N − 1]) où apparaı̂t chaque suffixe.

Critique : Comme le script demande de stocker en mémoire les séquences concaténées, la mémoire peut rapidement saturer
 lorsqu'on le fait tourner sur de plus grand groupe de séquences ou des séquences plus longues.
 Pistes d'améliorations a tester prochainement:
 - convertir en bytes ?
 - 

"""


# Import des packages nécessaires ---------------------------------------------------------------
from pydivsufsort import divsufsort
import sys
import numpy as np

# Lecture des séquences provenant du FOF (file of files) ----------------------------------------
def read_fof(fof_path):
    """
    Lit le FOF et retourne la liste des chemins de fichiers.
    
    Entrée : FOF 
    Sortie : liste des chemins de fichiers.

    """
    with open(fof_path) as f:
        return [line.strip() for line in f if line.strip()]


def read_fasta(path):
    """
    Lit un FASTA et retourne la séquence en majuscules, restreinte à ACGT.
    
    Entrée : Fichier FASTA indiqué par le FOF
    Sortie : Séquence contenu dans le fichier FASTA

    Notes : problème sur des grandes séquences, stockage en mémoire
    
    """
    parts = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                parts.append(line.strip().upper())
                
    return "".join(c for c in "".join(parts) if c in "ACGT") #sécurité pour s'assurer d'avoir bien une séquence


def build_from_fof(fof_path):
    """
    Concatène tous les génomes du FOF, et rajoute "$" à la fin de chacun d'entre eux.

    Entrée : Path du fof fournit.
    Sorties: 
    - (str) parts = tout les génomes concaténés du FOF, séparé par un "$".
    - (list) genome_of[i] = indice du génome auquel appartient la position i.
    - (list) paths = chemin des génomes, un indice de cette liste = la couleur (color arrray)
    """
    paths = read_fof(fof_path)
    parts = []
    genome_of = []

    for i, path in enumerate(paths):
        seq = read_fasta(path) + "$"
        parts.append(seq)
        genome_of.extend([i] * len(seq))

    return "".join(parts), genome_of, paths



# Réutilisation des fonctions du TP BWT ----------------------------------------------------------


### Construction transformée BWT 
def get_bwt(sequence : str, suffix_array : list[int]):
    """
    Construit la BWT à partir d'une séquence et de son suffix array.

    Pour chaque valeur SA[i], on prend le caractère qui précède le suffix dans la séquence

    Entrées :
    - sequence (str) finie par $
    - suffix array de la séquence

    Sortie: (str) Transformée de Burrows-Wheeler.
    """
    bwt = []

    for position in suffix_array:
        if position == 0:
            bwt.append("$")
        else:
            bwt.append(sequence[position - 1])

    return "".join(bwt)


### Construction du tableau C / N
def get_n(bwt):
    """
    Renvoi un dictionnaire avec la première occurance de la lettre dans la bwt et
    construit le tableau N utilisé dans le calcul LF.

    N[c] correspondant au nombre de caractères strictement inférieurs
    à c dans '$ACGT'.

    Entrée: (str) bwt
    Sortir : (dict) tableau N
    """
    nb_occ = {letter : 0 for letter in 'ATCG$'}
    for letter in bwt :
        nb_occ[letter]+=1

    n= {letter : 0 for letter in 'ATGC$'}
    n['$'] =0
    n['A'] =nb_occ['$']
    n['C']= n['A'] + nb_occ['A']
    n['G']= n['C'] + nb_occ['C']
    n['T']= n['G'] + nb_occ['G']
    return n


### Construction du tableau Rank
def get_r(bwt):
    """
    Construit le tableau Rank qui renvoi le rang sous forme de tableau
    r[i] = rang dans bwt du caractère bwt[i].

    Entrée : (str) bwt
    Sortie: (list[int]) tableau des rangs
    """

    dico = {'$':0, 'A':0, 'C':0, 'G':0, 'T':0}
    tab = list()
    for i,letter in enumerate(bwt):
        dico[letter] += 1
        tab.append(dico[letter])
        #print(f'{dico[letter]}\t{letter}')
            
    return  tab


### Fonction left_firt
def left_first(char, rank, n):
    """
    Calcule la position dans la colonne F correspondant à une occurrence
    donnée dans la colonne L (BWT).

    Entrées :
     - (str) caractère recherché
     - (int) rang du caractère
     - (dict) tableau obtenu avec get_n()

    Sortie : (int) possition du caractère dans la première colone

    """
    return n[char] + rank - 1


# Reconstruction de la séquence
def bwt_to_sequence(bwt, n, rank):
    """
    Reconstruit la séquence originale à partir de la BWT.

    Entrées: 
    - (str) bwt 
    - (dict) n provenant de get_n()
    - (list[int]) provenant de get_r()

    Sortie :
    - (str) séquence originale
    """
    sequence=""
    line =0
    while len(sequence) < len(bwt):
        sequence = bwt[line] + sequence
        line = left_first(bwt[line], rank[line], n)
    return sequence


def find_first(bwt, start, stop, letter):
    """Première ligne de [start, stop] dont la BWT vaut letter, sinon -1."""
    for line in range(start, stop + 1):
        if bwt[line] == letter:
            return line
    return -1


def find_last(bwt, start, stop, letter):
    """Dernière ligne de [start, stop] dont la BWT vaut letter, sinon -1."""
    for line in range(stop, start - 1, -1):
        if bwt[line] == letter:
            return line
    return -1


def backward_search(p, bwt, n, r):
    start, stop = 0, len(bwt) - 1

    for letter in reversed(p):
        first = find_first(bwt, start, stop, letter)
        if first == -1:
            return False, -1, -1
        last = find_last(bwt, start, stop, letter)

        start = left_first(letter, r[first], n)
        stop = left_first(letter, r[last], n)

    return True, start, stop


def get_occurrence(bwt, char, position):
    """
    Retourne le nombre d'occurrences de char dans BWT[:position].
    """
    return bwt[:position].count(char)


def contains(pattern, bwt, n, rank):
    """
    Indique si une séquence est présente.

    Sortie : Bolléen
    """
    found, _, _ = backward_search(
        pattern,
        bwt,
        n,
        rank
    )

    return found


def nb_occurrences(pattern, bwt, n, rank):
    """
    Retourne le nombre d'occurrences du motif.
    """
    found, start, stop = backward_search(
        pattern,
        bwt,
        n,
        rank
    )

    if not found:
        return 0

    return stop - start + 1




#---------------------------------------------------------------------------------------------
#                                    PROGRAMME PRINCIPAL
#---------------------------------------------------------------------------------------------

if __name__ == "__main__":

    # Demande le nombre d'arguments requis pour faire fonctionner le programme  
    if len(sys.argv) < 3:
        sys.exit("Usage : python script.py <fichier.fof> [motif]")

    # Initialisation
    fof_path = sys.argv[1]
    pattern = sys.argv[2]

    #---------------
    # Lecture des génomes et concaténation
    sequence, genome_of, paths = build_from_fof(fof_path)

    # Suffix array, BWT, N, Rank
    #Passage en bytes pour tenter de solutionner les pb de mémoires
    # seq_bytes = sequence.encode("ascii")
    # suffix_array = divsufsort(np.frombuffer(seq_bytes, dtype=np.uint8))
    suffix_array = divsufsort(sequence)
    bwt = get_bwt(sequence, suffix_array)
    n = get_n(bwt)
    rank = get_r(bwt)

    # Color array : génome d'origine de chaque suffixe
    color = [genome_of[pos] for pos in suffix_array]

    # Recherche du motif
    found, start, stop = backward_search(pattern, bwt, n, rank)

    #------------
    # Affichage / renvoit de programme
    print("Nb génomes: ", len(paths))
    print("Motif recherché: ", pattern)
    print("Présent: ", found)

    if found:
        genomes = sorted(set(color[start:stop + 1]))
        print("Occurrences: ", stop - start + 1)
        print("Génomes contenant le motif: ", [paths[g] for g in genomes])

