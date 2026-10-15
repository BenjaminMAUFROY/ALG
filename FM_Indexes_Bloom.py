# Modules
from pydivsufsort import divsufsort
import mmh3
from bitarray import bitarray
import os


# Classe FM_Index
class FM_Index:

    # Initialisation
    def __init__(self, seq):
        self._seq = seq
        self._sa = self.create_sa()
        self._bwt = self.create_bwt()
        self._n = self.create_n()
        self._r = self.create_r()
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
        return motif in self.bwt_2_seq()


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


# Fonction d'indexation de génome en N FM_Index et N filtres de Bloom
def indexation_2(path_or_filename: str, file_or_folder = "file"):
    if file_or_folder == "file":
        with open(path_or_filename) as filin:
            genome_list = filin.readlines().strip()
    elif file_or_folder == "folder":
        file_list = os.listdir(path_or_filename)
        genome_list = list()
        for file in file_list:
            with open(path_or_filename + file) as filin:
                genome_list.append()

    



if __name__ == "__main__":

    # Test classe FM_Index
    seq = FM_Index("GGCGGCACCGC$")
    sa = seq.get_sa()
    bwt = seq.get_bwt()
    r = seq.get_r()
    n = seq.get_n()
    fm_index = seq.fm_index()
    print(sa)
    print(bwt)
    print(r)
    print(n)
    print(fm_index)
    print(seq.contains("GGC"))

    # Test classe filtres de Bloom
    filtre = BloomFilter(7, 100000000)
    filtre.add("ACTTTGACC")
    for i in range(100000000):
        if filtre.get_matrice()[i]:
            print(f"my table {i} is true")
    filtre.add_kmers_from_sequence("acatactacgacgtacgatcgacgatcgatctatcagcgttggctactgatgctagtctgtttagcgactgtacg", 5)

