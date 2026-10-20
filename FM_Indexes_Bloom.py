# Modules
from pydivsufsort import divsufsort
import mmh3
from bitarray import bitarray
import math

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

    # Recherche un caractère dans la BWT en partant du debut
    def find_first(self, start, stop, letter):
        """Première ligne de [start, stop] dont la BWT vaut letter, sinon -1."""
        for line in range(start, stop + 1):
            if self.get_bwt()[line] == letter:
                return line
        return -1

    # Recherche un caractere dans la BWT en partant de la fin
    def find_last(self, start, stop, letter):
        """Dernière ligne de [start, stop] dont la BWT vaut letter, sinon -1."""
        for line in range(stop, start - 1, -1):
            if self.get_bwt()[line] == letter:
                return line
        return -1

    # Recherche un motif dans la seq
    def backward_search(self, motif):
        start = 0
        stop = len(self.get_bwt()) - 1

        for letter in reversed(motif):
            first = self.find_first(self.get_bwt(), start, stop, letter)
            if first == -1:
                return False
            last = self.find_last(self.get_bwt(), start, stop, letter)

            start = self.get_n()[self.get_bwt()[letter]] + self.get_r()[first] - 1
            stop = self.get_n()[self.get_bwt()[letter]] + self.get_r()[last] - 1

        return True

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
            with open(line) as seq:
                for line in seq:
                    if not line.startswith(">"):
                        genome += line.strip()
                n_kmer = len(genome)
                m = bf_size(n_kmer, h, fp_rate)
                filtre_bloom_genome = BloomFilter(h, m).add_kmers_from_sequence(genome, kmer)
                fm_index_genome = FM_Index(genome).fm_index()
                index.append({"Filtre_Bloom": filtre_bloom_genome, "FM_index": fm_index_genome,})
   
    return

    # Serialisation
    
        # Temps de serialisation
    
    # pkl.dump(index, open("genomes_indexees.dumped", "wb"))
 
if __name__ == "__main__":

    # Test classe FM_Index
    # seq = FM_Index("GGCGGCACCGC$")
    # sa = seq.get_sa()
    # bwt = seq.get_bwt()
    # r = seq.get_r()
    # n = seq.get_n()
    # fm_index = seq.fm_index()
    # print(sa)
    # print(bwt)
    # print(r)
    # print(n)
    # print(fm_index)
    # print(seq.contains("GGC"))

    # # Test classe filtres de Bloom
    # filtre = BloomFilter(7, 100000000)
    # filtre.add("ACTTTGACC")
    # for i in range(100000000):
    #     if filtre.get_matrice()[i]:
    #         print(f"my table {i} is true")
    # filtre.add_kmers_from_sequence("acatactacgacgtacgatcgacgatcgatctatcagcgttggctactgatgctagtctgtttagcgactgtacg", 5)

    # Test indexation_2
    # index = bf("./queries_coli_10k_len100bp.fa", 10000000)
    # print(len(index))
