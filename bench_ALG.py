"""Comparaison du temps de construction et de sérialisation
des méthodes generalized VS bf, avec un graphique barplot.

Utilisation : 
python3 bench_ALG.py genomes.fof

"""

#Modules
import re
import subprocess
import sys
import matplotlib.pyplot as plt

# ---- PARAMETRES ----------------------------------------------------------
SCRIPT = "index_colored_fm_index_V2.py"    # programme d'indexation testé
FOF = sys.argv[1]                          # fichier FoF à mentionné dans la commande
NB_REPETITIONS = 10                        # nombre de mesures par méthode
# --------------------------------------------------------------------------

METHODES = ["generalized", "bf"]


def mesurer(methode):
    """Lance l'indexation et renvoie (temps_build, temps_serialisation).

    Renvoie None si le programme échoue.
    """
    commande = [sys.executable, SCRIPT, "-i", FOF, "-m", methode, "-o", "bench_tmp.pkl"]
    resultat = subprocess.run(commande, capture_output=True, text=True)

    # On cherche les deux lignes affichées par le programme d'indexation
    build = re.search(r"OUT TIME_BUILD:\s*([0-9.eE+-]+)", resultat.stdout)
    serial = re.search(r"OUT TIME_SERIALISATION:\s*([0-9.eE+-]+)", resultat.stdout)

    if resultat.returncode != 0 or build is None or serial is None:
        print(f"Échec pour la méthode {methode} :")
        print(resultat.stderr.strip().splitlines()[-1:])
        return None

    return float(build.group(1)), float(serial.group(1))


# Mesures : répète NB_REPETITIONS fois et garde la moyenne
temps_build = {}
temps_serial = {}

for methode in METHODES:
    mesures = []
    for _ in range(NB_REPETITIONS):
        mesure = mesurer(methode)
        if mesure is None:
            mesures = []
            break   # inutile de répéter un échec
        mesures.append(mesure)

    if mesures:
        temps_build[methode] = sum(m[0] for m in mesures) / len(mesures)
        temps_serial[methode] = sum(m[1] for m in mesures) / len(mesures)
        print(f"{methode} : build {temps_build[methode]:.4f} s, "
              f"sérialisation {temps_serial[methode]:.4f} s")

if not temps_build:
    sys.exit("Aucune mesure réussie : rien à tracer.")


# -------- VISUALISATION ----------------
# Graphique : deux petits graphiques côte à côte
figure, (axe1, axe2) = plt.subplots(1, 2, figsize=(9, 4))

axe1.bar(list(temps_build.keys()), list(temps_build.values()))
axe1.set_title("Temps de construction")
axe1.set_ylabel("Temps (s)")

axe2.bar(list(temps_serial.keys()), list(temps_serial.values()))
axe2.set_title("Temps de sérialisation")
axe2.set_ylabel("Temps (s)")

figure.tight_layout()
figure.savefig("comparaison.png")
print("Graphique enregistré : comparaison.png")