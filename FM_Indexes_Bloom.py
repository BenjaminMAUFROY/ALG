import mmh3
from bitarray import bitarray

class BloomFilter:
    def __init__(self, H, M):
        self._H = H
        self._M = M
        matrice = bitarray(M)
        matrice.setall(0)
        self._matrice = matrice


    def get_H(self):
        return self._H
    def get_matrice(self):
        return self._matrice
    def get_M(self):
        return self._M

    def __contains__(self, item):
        return item in self.get_matrice()
    
    def add(self, item):
        for seed_id in range(self.get_H()):
            key = mmh3.hash(item, seed_id, signed=False) % self.get_M()
            self.get_matrice()[key] = 1

    def add_kmers_from_sequence(self, seq, k):
        for i in range(len(seq)-k+1):
            self.add(seq[i:i+k])
        return self.get_matrice()
        


if __name__ == "__main__":
    filtre = BloomFilter(7, 100000000)
    filtre.add("ACTTTGACC")
    # for i in range(100000000):
    #     if filtre.get_matrice()[i]:
    #         print(f"my table {i} is true")
    filtre.add_kmers_from_sequence("acatactacgacgtacgatcgacgatcgatctatcagcgttggctactgatgctagtctgtttagcgactgtacg", 5)