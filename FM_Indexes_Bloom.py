"""import mmh3
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
        
class 



if __name__ == "__main__":
    filtre = BloomFilter(7, 100000000)
    filtre.add("ACTTTGACC")
    # for i in range(100000000):
    #     if filtre.get_matrice()[i]:
    #         print(f"my table {i} is true")
    filtre.add_kmers_from_sequence("acatactacgacgtacgatcgacgatcgatctatcagcgttggctactgatgctagtctgtttagcgactgtacg", 5)
"""

"""
______________________________________FM-Index____________________________________________________________
"""
from pydivsufsort import divsufsort

b = "GGCGGCACCGC$"
ba = divsufsort(b)
print("i\tSA\tf\tSuffixes")
for i in range(len(b)):
    print(f"{i}\t{ba[i]}\t{b[ba[i]]}\t{b[ba[i]:]}")

def get_bwt(b,ba):
    return [b[i-1] for i in ba]

def get_n(bwt):
    count = {"A":0, "C":0, "G":0, "$":0}
    n = {"$":0}
    for letter in bwt:
        count[letter] += 1

    n["A"] = count["$"]
    n["C"] = count["$"] + count["A"]
    n["G"] = count["$"] + count["A"] + count["C"]

    n["T"] = count["$"] + count["A"] + count["C"] + count["G"]

    return n

def get_r(bwt):
    count = {"A":0, "C":0, "G":0, "T":0, "$":0}
    r = []
    for letter in bwt:
        count[letter] += 1
        r.append(count[letter])
    return r

def left_first(a, k, n):
    return n[a] + k - 1

def bwt_2_seq(bwt, n, r):
    seq = ["$",]
    count = 0
    for _ in range(len(bwt)-1):
        seq.append(bwt[count])
        count = left_first(bwt[count], r[count], n)
    return "".join(seq[::-1])

def contains(p, bwt, n, r):

    return p in bwt_2_seq(bwt, n, r)



print(get_n(get_bwt(b,ba)))
print(get_r(get_bwt(b, ba)))
bwt = get_bwt(b, ba)
r = get_r(bwt)
n = get_n(bwt)
print(bwt_2_seq(bwt, n, r))

# le FM-index contient l'ensemble des suffixes (f), l'index de l'ordre alphabetique des suffixes (i) et la BWT, sachant qu'on ne stock pas l'entièreté des suffixes