from UniprotID_GO import uniprotid_to_go
from ast import literal_eval
import collections


def func(line):
    ts = line.split("\t")
    tsv = ts[298][0:-1]
    k = str(int(ts[0].replace("OG", "")))
    if tsv != "":
        ids = [uniprotid_to_go(header.split("_")[1])
               for header in tsv.split(", ")]
        ids = [a for b in ids for a in b]
        ids_c = [id for id in ids if "C:" in id]
        ids_f = [id for id in ids if "F:" in id]
        ids_p = [id for id in ids if "P:" in id]
        ids_c = collections.Counter(ids_c).most_common()
        ids_f = collections.Counter(ids_f).most_common()
        ids_p = collections.Counter(ids_p).most_common()
        ids_c = ','.join(map(str, ids_c))
        ids_f = ','.join(map(str, ids_f))
        ids_p = ','.join(map(str, ids_p))
        res = "\t".join([str(k), ids_c, ids_f, ids_p])

        return res
    else:
        return "\t".join([str(k), "", "", ""])


