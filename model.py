# Modules
from pydivsufsort import divsufsort
import mmh3
from bitarray import bitarray
import math

# Classe FM_Index
class FM_Index:

    # Initialisation
    def __init__(self, seq, genome_of=None):
        self._seq = seq
        self._sa = self.create_sa()
        self._bwt = self.create_bwt()
        self._n = self.create_n()
        self._r = self.create_r()
        self._color = self.create_color(genome_of)
        self._fm_index = self.fm_index()

    # Méthodes
    # Création du suffixe_array
    def create_sa(self):
        return divsufsort(self.get_seq())

    # Création de la BWT
    def create_bwt(self):
        return [self.get_seq()[i-1] for i in self.get_sa()]

    # Création de n
    def create_n(self):
        count = {"A":0, "C":0, "G":0, "T":0, "$":0,}
        n = {"$":0}
        for letter in self.get_bwt():
            count[letter] += 1

        n["A"] = count["$"]
        n["C"] = count["$"] + count["A"]
        n["G"] = count["$"] + count["A"] + count["C"]
        n["T"] = count["$"] + count["A"] + count["C"] + count["G"]

        return n

    # Création de r
    def create_r(self):
        count = {"A":0, "C":0, "G":0, "T":0, "$":0}
        r = []
        for letter in self.get_bwt():
            count[letter] += 1
            r.append(count[letter])
        return r

    def create_color(self, genome_of):
        """
        Construit le color array : génome d'origine de chaque suffixe, aligné sur le suffix array.
        
        Entrée : liste genome_of, avec la couleur de chaque position de la séquence
        (None pour un FM_Index sans color array)
        Sortie : np.ndarray color array 
        (ou None si genome_of est None)
        """
        if genome_of is None:
            return None
        return [genome_of[position] for position in self.get_sa()]
    
    # FM_Index
    def fm_index(self):
        return {"bwt": self.get_bwt(), "n": self.get_n(), "r": self.get_r(),}

    # Reconstruit la sequence à partir de la BWT
    def bwt_2_seq(self):
        seq = ["$",]
        count = 0
        for _ in range(len(self.get_bwt())-1):
            seq.append(self.get_bwt()[count])
            count = self.get_n()[self.get_bwt()[count]] + self.get_r()[count] - 1
        return "".join(seq[::-1])

    
    # Compare un motif à la sequence indexée
    def contains(self, motif):
        """
        Indique si une séquence est présente.

        Sortie : Bolléen
        """
        return self.backward_search(motif)


    # Définition des accesseurs
    def get_seq(self):
        return self._seq
    def get_bwt(self):
        return self._bwt
    def get_n(self):
        return self._n
    def get_r(self):
        return self._r
    def get_sa(self):
        return self._sa
    def get_fm_index(self):
        return self._fm_index
    def get_color(self):
        return self._color


class BloomFilter:

    # Initialisation
    def __init__(self, H, M):
        self._H = H
        self._M = M
        matrice = bitarray(M)
        matrice.setall(0)
        self._matrice = matrice

    # Méthodes
    # Surcharge de la classe contains pour pouvoir effectuer des tests de comparaison
    def __contains__(self, item):
        return item in self.get_matrice()

    # Création de add qui remplit la matrice de bit de 1 en fonction de la sortie des fonctions de hashage sur un kmer item    
    def add(self, item):
        for seed_id in range(self.get_H()):
            key = mmh3.hash(item, seed_id, signed=False) % self.get_M()
            self.get_matrice()[key] = 1

    # Méthode qui renvoie une matrice de Bloom remplie à partir d'une séquence qui est découpée en kmer de taille k  
    def add_kmers_from_sequence(self, seq, k):
        for i in range(len(seq)-k+1):
            self.add(seq[i:i+k])
        return self.get_matrice()

    # Définition des accesseurs
    def get_H(self):
        return self._H
    def get_matrice(self):
        return self._matrice
    def get_M(self):
        return self._M

# Fonction pour la lecture du FOF
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

    Entrée : Adresse du fof fournit.
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

# Fonction de calcul de la taille de la matrice du filtre de bloom
def bf_size(n_kmer: int, H = 7, fp_rate = 0.01) -> int:
    return  math.ceil(- H * (n_kmer / (math.log(1 - fp_rate ** (1 / H)))))

# Fonction d'indexation de génome en N FM_Index et N filtres de Bloom
def bf(FOF_name: str) -> list:

    # Parametres des filtres de bloom et fm_index
    h = 7           # Nombre de fonction de hashages
    fp_rate = 0.01  # Ratio de faux positifs tolérés
    kmer = 31       # Taille des kmers

    # Indexation
    index = []
    with open(FOF_name) as filin:
        for line in filin:
            genome = ""
            with open(line.strip()) as seq:
                for line in seq:
                    if not line.startswith(">"):
                        genome += line.strip()
                n_kmer = len(genome)
                m = bf_size(n_kmer, h, fp_rate)
                filtre_bloom_genome = BloomFilter(h, m).add_kmers_from_sequence(genome, kmer)
                fm_index_genome = FM_Index(genome).fm_index()
                index.append({"Filtre_Bloom": filtre_bloom_genome, "FM_index": fm_index_genome,})
   
    return index

# Méthode de requêtage d'un motif dans un filtre de Bloom
def query(matrice, item, h = 7):
    for seed_id in range(h):
        key = mmh3.hash(item, seed_id, signed=False) % len(matrice)
        if matrice[key]:
            return True
    return False

# Recherche un caractère dans la BWT en partant du debut
def find_first(bwt, start, stop, letter):
    """Première ligne de [start, stop] dont la BWT vaut letter, sinon -1."""
    for line in range(start, stop + 1):
        if bwt[int(line)] == letter:
            return line
    return -1

# Recherche un caractere dans la BWT en partant de la fin
def find_last(bwt, start, stop, letter):
    """Dernière ligne de [start, stop] dont la BWT vaut letter, sinon -1."""
    for line in range(stop, start - 1, -1):
        if bwt[int(line)] == letter:
            return line
    return -1

# Recherche un motif dans la seq
def backward_search(fm_index, motif):
    start = 0
    stop = len(fm_index["bwt"]) - 1

    for letter in reversed(motif):
        first = find_first(fm_index["bwt"], start, stop, letter)
        if first == -1:
            return False
        last = find_last(fm_index["bwt"], start, stop, letter)

        start = fm_index["n"][letter] + fm_index["r"][first] - 1
        stop = fm_index["n"][letter] + fm_index["r"][last] - 1

    return True
