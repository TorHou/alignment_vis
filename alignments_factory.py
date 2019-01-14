import subprocess
from argparse import ArgumentParser
from Bio import SeqIO
import matplotlib.pyplot as plt
import pickle

parser = ArgumentParser()
parser.add_argument("contigfile", help="Contig File")
parser.add_argument("contigpairs", help="Overlapping Contig-Pairs File")

args = parser.parse_args()

pairs = []
with open(args.contigpairs) as f:
    for line in f:
        if not line.startswith(">"):
            c1, c2 = line.rstrip().split("\t")
            pairs.append((c1,c2))

contigs = {}
for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = str(read.seq)

#pairs = ["1017APD_1409APD", "1036APD_1088APD", "1045APD_1377APD", "1727APD_1110APD", "1144APD_1686APD", "1176APD_1733APD"]
#pairs = ["1017APD_1409APD"]
#pairs.append
#pairs.append("613APD_1758APD")
#pairs.append("691APD_1790APD")

pair_scores = {}

def shortname(ctg):
    if "_" in ctg:
        return ctg.split("_")[1]
    else:
        return ctg

for i, pair in enumerate(pairs):
    print(str(i) + ": " + str(pair))
    ctg1, ctg2 = pair
    with open("tmp1.fasta", "w") as f:
        f.write(">" + ctg1 +"\n")
        f.write(contigs[shortname(ctg1)] + "\n")
    with open("tmp2.fasta", "w") as f:
        f.write(">" + ctg2 +"\n")
        f.write(contigs[shortname(ctg2)] + "\n")
    points = subprocess.run(["./alignments", "tmp1.fasta", "tmp2.fasta"], stdout=subprocess.PIPE, universal_newlines=True)
    pair_scores[pair] = points.stdout.rstrip().split("\n")
    
    #print(points["stdout"])
    p2 = [float(x) for x in points.stdout.rstrip().split("\n")]
    maxov = len(p2)
    plt.figure(i)
    plt.axis([-200, 0, -1.05, 1.05])
    plt.plot(range(-maxov,0), p2)
    plt.xlabel(str(pair))
    plt.ylabel("normalized score")
    plt.savefig(shortname(ctg1) + "_" + shortname(ctg2) + ".png")
    plt.close()
    
    
pickle.dump( pair_scores, open("pair_scores.pkl", "wb"))
