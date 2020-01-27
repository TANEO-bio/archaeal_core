import sys
import math
from UniProt_to_GO import term
from UniprotID_GO import uniprotid_to_GO_Term


def func(line):
    ts = line.split("\t")
    og = ts[0]
    print(len(ts))
    tsv = ts[298]
    if tsv != "":
        ids = {header.split("_")[1] for header in tsv.split(", ")}
        tmp = [term(uniprotid_to_GO_Term(i)[0]).tree() if len(
            uniprotid_to_GO_Term(i)) > 0 else [] for i in ids]
        tmp = [a for b in tmp for a in b]
        tmp = set(tmp)
        return str(int(og.replace("OG", ""))) + "\t" + og + "\t" + "@".join(tmp)
    else:
        return str(int(og.replace("OG", ""))) + "\t" + og + "\t" + ""


num = 0
f = open("Orthogroups_test.tsv", "r")
lines = f.readlines()[1:]
f.close()

lines = lines[num::10]

for line in lines:
    print(line)
    print(str(func(line)))
