import glob
import random
import subprocess
from ete3 import NCBITaxa


ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def taxid2lineage(taxid):
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    lineage_dict = dict()
    names = ncbi.get_taxid_translator(lineage)
    for rank in ranks:
        for k, v in ncbi.get_rank(lineage).items():
            if v == rank:
                lineage_dict.update({v: names[k]})
    return lineage_dict


path_list = glob.glob("./data/*.fasta")
path_list.remove('./data/swiss_prot_archaea_2.fasta')

print("dict_tax = {", end="")

for path in path_list:
    cmd = "cat {} | head -n 1 | awk -F '@' '{{print $3}}'".format(path)
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    taxid = proc.stdout.decode("utf8").replace("\n", "")
    taxonomy_dict = taxid2lineage(taxid)
    taxonomy = []
    for rank in ranks:
        try:
            taxonomy.append(taxonomy_dict[rank])
        except KeyError:
            taxonomy.append("unknown{}".format(random.randint(1000000000, 9999999999)))
    taxonomy = "/".join(taxonomy)
    taxonomy = taxonomy.replace(" ", "_")
    print("\"" + str(taxid) + "\": " + "\"" + taxonomy + "/\", ", end="")
print("}")
