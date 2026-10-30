"""Modèle d'indexation de génomes : FM-index et filtres de Bloom.

Ce module fournit les classes et fonctions utilisées par le script
d'indexation :

- ``FM_Index`` : construction d'un FM-index (suffix array, BWT, tables n et r)
  à partir d'une séquence terminée par "$", avec un color array optionnel qui
  associe à chaque suffixe le génome d'origine.
- ``BloomFilter`` : filtre de Bloom de k-mers reposant sur la fonction de
  hachage mmh3 et sur un tableau de bits (bitarray).
- Des fonctions de lecture (FOF, FASTA), de dimensionnement des filtres de
  Bloom, d'indexation (``bf``) et de requêtage (``query``, ``backward_search``,
  ``backward_search_interval``).

Deux stratégies d'indexation sont visées :

- ``build_from_fof`` : concaténation de tous les génomes (chacun terminé par
  "$") en une seule séquence, destinée à un FM-index unique accompagné d'un
  color array ;
- ``bf`` : un FM-index et un filtre de Bloom par génome.

Les paramètres du filtre de Bloom (nombre de fonctions de hachage h, taux de
faux positifs fp_rate, taille des k-mers kmer) sont des arguments de ``bf``.
Au requêtage, ``query`` doit être appelée avec le même nombre de fonctions de
hachage h, sur des k-mers de la même taille que ceux insérés dans le filtre.
"""

# Modules
import math

from pydivsufsort import divsufsort
import mmh3
from bitarray import bitarray


# Classe FM_Index
class FM_Index:
    """FM-index d'une séquence nucléotidique.

    La séquence doit être composée des lettres A, C, G et T et se terminer par
    le caractère "$". Plusieurs "$" sont possibles dans le cas d'une
    concaténation de génomes.

    Attributs :
    - _seq (str) : séquence indexée.
    - _sa : suffix array de la séquence (positions de départ des suffixes
      triés).
    - _bwt (list) : transformée de Burrows-Wheeler, liste de caractères.
    - _n (dict) : pour chaque lettre, nombre de caractères de la séquence
      strictement inférieurs à cette lettre.
    - _r (list) : pour chaque position de la BWT, rang (à partir de 1) de son
      caractère parmi les occurrences identiques jusqu'à cette position.
    - _color (list ou None) : color array, aligné sur le suffix array.
    - _fm_index (dict) : regroupe la BWT et les tables n et r.
    """

    def __init__(self, seq, genome_of=None):
        """
        Construit le suffix array, la BWT, les tables n et r, le color array
        et le dictionnaire FM-index.

        Entrée :
        - (str) seq : séquence à indexer, terminée par "$".
        - (list) genome_of : génome d'origine de chaque position de seq
          (None pour un FM_Index sans color array).
        """
        self._seq = seq
        self._sa = self.create_sa()
        self._bwt = self.create_bwt()
        # Tables n et r, calculées à partir de la BWT
        self._n = self.create_n()
        self._r = self.create_r()
        # Color array (None si genome_of n'est pas fourni)
        self._color = self.create_color(genome_of)
        # Regroupement de la BWT et des tables n et r dans un dictionnaire
        self._fm_index = self.fm_index()

    # Méthodes
    def create_sa(self):
        """
        Calcule le suffix array de la séquence avec divsufsort.

        Sortie : tableau numpy des positions de départ des suffixes,
        triés par ordre lexicographique.
        """
        return divsufsort(self.get_seq())

    def create_bwt(self):
        """
        Construit la transformée de Burrows-Wheeler (BWT).

        Pour chaque suffixe du suffix array, on prend le caractère qui le
        précède dans la séquence. Pour le suffixe de position 0, l'indice -1
        donne le dernier caractère de la séquence, c'est-à-dire "$".

        Sortie : (list) liste de caractères.
        """
        return [self.get_seq()[i - 1] for i in self.get_sa()]

    def create_n(self):
        """
        Construit la table n : pour chaque lettre, nombre de caractères de la
        séquence strictement inférieurs à cette lettre (ordre $ < A < C < G
        < T). n[lettre] est donc l'indice de la première ligne dont la
        première colonne de la matrice triée commence par cette lettre.

        Sortie : (dict) {"$": 0, "A": ..., "C": ..., "G": ..., "T": ...}
        """
        count = {"A": 0, "C": 0, "G": 0, "T": 0, "$": 0}
        n = {"$": 0}
        for letter in self.get_bwt():
            count[letter] += 1

        # Somme cumulée des effectifs des lettres plus petites
        n["A"] = count["$"]
        n["C"] = count["$"] + count["A"]
        n["G"] = count["$"] + count["A"] + count["C"]
        n["T"] = count["$"] + count["A"] + count["C"] + count["G"]

        return n

    def create_r(self):
        """
        Construit la table r : pour chaque position de la BWT, rang (à partir
        de 1) de son caractère parmi les occurrences de cette même lettre dans
        la BWT jusqu'à cette position incluse.

        Avec n, elle permet le LF-mapping : n[bwt[i]] + r[i] - 1.

        Sortie : (list) liste d'entiers de même longueur que la BWT.
        """
        count = {"A": 0, "C": 0, "G": 0, "T": 0, "$": 0}
        r = []
        for letter in self.get_bwt():
            count[letter] += 1
            r.append(count[letter])
        return r

    def create_color(self, genome_of):
        """
        Construit le color array : génome d'origine de chaque suffixe, aligné
        sur le suffix array.

        Entrée : liste genome_of, avec la couleur de chaque position de la
        séquence (None pour un FM_Index sans color array)
        Sortie : (list) color array, c'est-à-dire genome_of[position] pour
        chaque position du suffix array (ou None si genome_of est None)
        """
        # Pas de color array demandé
        if genome_of is None:
            return None

        # Couleur de chaque suffixe = couleur de sa position de départ
        return [genome_of[position] for position in self.get_sa()]

    def fm_index(self):
        """
        Regroupe dans un FM-index les éléments nécessaires à la recherche
        de motifs.

        Sortie : (dict) {"bwt": BWT, "n": table n, "r": table r}.
        """
        return {"bwt": self.get_bwt(), "n": self.get_n(), "r": self.get_r()}

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
    """Filtre de Bloom destiné à stocker les k-mers d'une séquence.

    Un élément ajouté est haché H fois par mmh3 et les H positions obtenues,
    modulo M, sont mises à 1 dans un tableau de bits.

    Attributs :
    - _H (int) : nombre de fonctions de hachage.
    - _M (int) : taille du tableau de bits.
    - _matrice (bitarray) : tableau de bits, initialisé à 0.
    """

    # Initialisation
    def __init__(self, H, M):
        """
        Créer un filtre de Bloom vide.

        Entrée :
        - (int) H : nombre de fonctions de hachage.
        - (int) M : taille du tableau de bits.
        """
        self._H = H
        self._M = M
        # Tableau de bits de taille M, entièrement mis à 0
        matrice = bitarray(M)
        matrice.setall(0)
        self._matrice = matrice

    # Méthodes
    def add(self, item):
        """
        Ajoute un élément (un k-mer) au filtre.

        Pour chacune des H fonctions, l'élément est haché puis la position
        obtenue (modulo M) est mise à 1 dans le tableau de bits.

        Entrée : (str) item, k-mer à ajouter.
        """
        for seed_id in range(self.get_H()):
            # Position dans le tableau de bits pour cette graine
            key = mmh3.hash(item, seed_id, signed=False) % self.get_M()
            self.get_matrice()[key] = 1

    def add_kmers_from_sequence(self, seq, k):
        """
        Découpe une séquence en k-mers (fenêtre glissante de pas 1) et les
        ajoute tous au filtre.

        Entrée :
        - (str) seq : séquence à découper.
        - (int) k : taille des k-mers.
        Sortie : (bitarray) tableau de bits du filtre, une fois rempli
        (et non l'objet BloomFilter lui-même).
        """
        for i in range(len(seq) - k + 1):
            self.add(seq[i:i + k])
        return self.get_matrice()

    # Définition des accesseurs
    def get_H(self):
        return self._H

    def get_matrice(self):
        return self._matrice

    def get_M(self):
        return self._M


def read_fof(fof_path):
    """
    Lit le FOF et retourne la liste des chemins de fichiers.

    Entrée : FOF
    Sortie : liste des chemins de fichiers.

    """
    # Une ligne = un chemin ; les lignes vides sont ignorées
    with open(fof_path) as f:
        return [line.strip() for line in f if line.strip()]


def read_fasta(path):
    """
    Lit un FASTA et retourne la séquence en majuscules, restreinte à ACGT.

    Entrée : Fichier FASTA indiqué par le FOF.
    Sortie : Séquence contenue dans le fichier FASTA.

    Les lignes d'en-tête (">") sont ignorées : si le fichier contient
    plusieurs entrées, leurs séquences sont concaténées en une seule. Les
    caractères autres que A, C, G et T (N, par exemple) sont supprimés.

    """
    parts = []
    with open(path) as f:
        for line in f:
            # Les lignes d'en-tête (">") sont ignorées
            if not line.startswith(">"):
                parts.append(line.strip().upper())

    # Sécurité pour s'assurer d'avoir bien une séquence
    return "".join(c for c in "".join(parts) if c in "ACGT")


def build_from_fof(fof_path):
    """
    Concatène tous les génomes du FOF, et rajoute "$" à la fin de chacun
    d'entre eux.

    Entrée : Adresse du FOF fourni.
    Sorties:
    - (str) parts = tous les génomes concaténés du FOF, chacun terminé par un
      "$".
    - (list) genome_of[i] = indice du génome auquel appartient la position i.
    - (list) paths = chemins des génomes, l'indice dans cette liste est la
      couleur (color array).
    """
    paths = read_fof(fof_path)
    parts = []
    genome_of = []

    for i, path in enumerate(paths):
        seq = read_fasta(path) + "$"
        parts.append(seq)
        # Chaque position de ce génome (y compris le "$") reçoit la couleur i
        genome_of.extend([i] * len(seq))

    return "".join(parts), genome_of, paths


def bf_size(n_kmer: int, H=7, fp_rate=0.01) -> int:
    """
    Calcule la taille (en bits) d'un filtre de Bloom pour un taux de faux
    positifs donné.

    Formule : M = ceil(-H * n / ln(1 - p^(1/H))), avec n le nombre d'éléments
    à insérer, H le nombre de fonctions de hachage et p le taux de faux
    positifs.

    Entrée :
    - (int) n_kmer : nombre d'éléments à insérer.
    - (int) H : nombre de fonctions de hachage (7 par défaut).
    - (float) fp_rate : taux de faux positifs toléré (0.01 par défaut).
    Sortie : (int) taille du tableau de bits.
    """
    return math.ceil(-H * (n_kmer / (math.log(1 - fp_rate ** (1 / H)))))


def bf(FOF_name: str, h=7, fp_rate=0.01, kmer=31) -> list:
    """
    Indexe chaque génome du FOF sous la forme d'un filtre de Bloom des k-mers
    dudit génome et d'un FM-index.

    Les chemins du FOF sont lus avec read_fof (les lignes vides sont
    ignorées). Chaque génome est lu avec read_fasta puis terminé par un
    unique "$".

    Entrée :
    - (str) FOF_name : chemin du FOF (un chemin de fichier FASTA par ligne).
    - (int) h : nombre de fonctions de hachage du filtre de Bloom (7 par
      défaut).
    - (float) fp_rate : taux de faux positifs visé, en fraction (0.01 par
      défaut, soit 1 %) ; il fixe la taille du filtre via bf_size.
    - (int) kmer : taille des kmers insérés dans le filtre (31 par défaut).
    Au requêtage, query doit utiliser le même h et des kmers de la même
    taille.
    Sortie : (list) une entrée par génome, dans l'ordre du FOF :
    {"Filtre_Bloom": bitarray du filtre,
    "FM_index": dictionnaire {"bwt", "n", "r"}}.
    """

    # Indexation
    index = []
    # Une ligne du FOF = un génome
    for path in read_fof(FOF_name):
        # Lecture du génome : majuscules et restreint à ACGT. Un "$" est
        # ajouté à la fin du génome
        genome = read_fasta(path) + "$"
        # Dimensionnement du filtre de Bloom pour ce génome : n_kmer est le
        # nombre de sous-séquences du génome de taille kmer.
        n_kmer = len(genome) - kmer + 1
        m = bf_size(n_kmer, h, fp_rate)
        # Filtre de Bloom rempli avec les kmers du génome
        filtre_bloom_genome = BloomFilter(h, m).add_kmers_from_sequence(
            genome, kmer)
        # FM-index du génome
        fm_index_genome = FM_Index(genome).fm_index()
        index.append({"Filtre_Bloom": filtre_bloom_genome,
                      "FM_index": fm_index_genome})

    return index


def query(matrice, item, h=7):
    """
    Interroge un filtre de Bloom pour savoir si un kmer est présent.

    Le kmer est déclaré présent si les h positions données par les h
    fonctions de hachage valent toutes 1.

    Entrée :
    - (bitarray) matrice : tableau de bits du filtre.
    - (str) item : kmer recherché, de la même taille que les kmers insérés.
    - (int) h : nombre de fonctions de hachage, identique à celui utilisé à
      la construction.
    Sortie : Booléen, False si le kmer est absent, True s'il est présent.
    """
    # Le kmer est absent dès que l'une des positions vaut 0
    for seed_id in range(h):
        key = mmh3.hash(item, seed_id, signed=False) % len(matrice)
        if not matrice[key]:
            return False
    return True


# 
def find_first(bwt, start, stop, letter):
    """
    Recherche un caractère dans la BWT en partant du debut. Première ligne de
    [start, stop] dont la BWT vaut letter, sinon -1.

    Entrée :
    - bwt (list) : transformée de Burrows-Wheeler, liste de caractères.
    - start (int) : position de début de l'intervalle.
    - stop (int) : position de fin de l'intervalle.
    - letter (str) : caractère recherché.
    Sortie : int, la position du caractère recherché.
    """
    for line in range(start, stop + 1):
        if bwt[int(line)] == letter:
            return line
    return -1


# 
def find_last(bwt, start, stop, letter):
    """
    Recherche un caractere dans la BWT en partant de la fin. Dernière ligne
    de [start, stop] dont la BWT vaut letter, sinon -1.

    Entrée :
    - bwt (list) : transformée de Burrows-Wheeler, liste de caractères.
    - start (int) : position de début de l'intervalle.
    - stop (int) : position de fin de l'intervalle.
    - letter (str) : caractère recherché.
    Sortie : int, la position du caractère recherché.
    """
    for line in range(stop, start - 1, -1):
        if bwt[int(line)] == letter:
            return line
    return -1


def backward_search(fm_index, motif):
    """
    Recherche un motif dans une séquence indexée par recherche arrière
    (backward search).

    Entrée :
    - (dict) fm_index : dictionnaire {"bwt", "n", "r"} (voir
      FM_Index.fm_index).
    - (str) motif : séquence à rechercher.
    Sortie : Booléen, True si le motif est présent dans la séquence indexée.
    """
    # Au départ, l'intervalle couvre toutes les lignes
    start = 0
    stop = len(fm_index["bwt"]) - 1

    # Lecture du motif de droite à gauche
    for letter in reversed(motif):
        first = find_first(fm_index["bwt"], start, stop, letter)
        if first == -1:
            # Si la lettre est absente de l'intervalle le motif n'est pas
            # présent
            return False
        last = find_last(fm_index["bwt"], start, stop, letter)

        # Nouvel intervalle par LF-mapping : n[lettre] + rang - 1
        start = fm_index["n"][letter] + fm_index["r"][first] - 1
        stop = fm_index["n"][letter] + fm_index["r"][last] - 1

    return True


def backward_search_interval(fm_index, motif):
    """
    Recherche un motif par recherche arrière (backward search) et retourne
    l'intervalle de lignes de la matrice triée des suffixes dont le début est
    le motif.

    Même parcours que backward_search, mais l'intervalle est retourné au lieu
    d'un booléen. Ses lignes sont aussi celles du suffix array et du color
    array : il permet de retrouver les génomes d'origine des occurrences
    (méthode generalized).

    Entrée :
    - (dict) fm_index : dictionnaire {"bwt", "n", "r"} (voir
      FM_Index.fm_index).
    - (str) motif : séquence à rechercher.
    Sortie : (tuple ou None) (start, stop), bornes incluses de l'intervalle de
    lignes où le motif est présent ; None si le motif est absent (ou s'il
    contient une lettre absente de l'index).
    """
    # Au départ, l'intervalle couvre toutes les lignes
    start = 0
    stop = len(fm_index["bwt"]) - 1

    # Lecture du motif de droite à gauche
    for letter in reversed(motif):
        # Lettre inconnue de l'index : le motif ne peut pas être présent
        if letter not in fm_index["n"]:
            return None

        first = find_first(fm_index["bwt"], start, stop, letter)
        if first == -1:
            return None

        last = find_last(fm_index["bwt"], start, stop, letter)

        # Nouvel intervalle par LF-mapping : n[lettre] + rang - 1
        start = fm_index["n"][letter] + fm_index["r"][first] - 1
        stop = fm_index["n"][letter] + fm_index["r"][last] - 1

        if start > stop:
            return None

    return start, stop
