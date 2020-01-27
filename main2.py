import subprocess
import pandas as pd
import collections
import scipy.stats as stats
# from UniprotID_GO import uniprotid_to_go
from Taxid_Taxonomy import taxid_to_taxonomy, get_taxonomy_list
from UniProt_GO import funcs
from dedup_taxa import taxon


f = open("Orthogroups.tsv")
lines = f.readlines()
f.close()

g = open("Taxid_Taxonomy.py")
txt = g.read()
g.close()

h = open("result.txt")
h_lines = h.readlines()
h.close()

ar_csv = open("ar14.arCOG.csv")
ar_csv_lines = ar_csv.readlines()
ar_csv.close()

ranks = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
len_tsv = sum([1 for _ in open('Orthogroups.tsv')]) - 1


class OrthoGroup():
    def __init__(self, number):
        og = lines[number + 1].split("\t")
        og = [header.replace("\n", "")
              for header in og if header != "" and "_" in header]
        self.OrthologGroup = og
        self.Number = number
        self.Headers = og[-1].replace(" ", "").split(",")
        if "@" not in self.Headers[0]:
            self.Name = "_".join(self.Headers[0].split("_")[
                                 4:]).split("_OS=")[0]
#            self.UniProtID = [header.split("_")[1] for header in self.Headers]
#            self.COGFunc = [
#                h_line.split("\t")[1].split("@")[1] for id_ in self.UniProtID for h_line in h_lines if id_ in h_line]
#            self.COGFunc = collections.Counter([ar_csv_line.split(
#                ",")[8] for id_ in self.COGFunc for ar_csv_line in ar_csv_lines if id_ in ar_csv_line if len(ar_csv_line.split(",")) > 8])
#            self.COGFunc = self.COGFunc.most_common()
        else:
            self.Name = og[-1].split("@")[5]
#            self.UniProtID = None

    def get_func(self):
        return num_to_func(self.Number)

    def get_conservation(self, taxonomy):
        headers = self.OrthologGroup
        taxids = [header.split("@")[2] for header in headers if "@" in header]
        taxonomies = [taxid_to_taxonomy(taxid) for taxid in taxids]
        num_conservation = len(
            [tax for tax in taxonomies if taxonomy + "/" in tax])
        num_all_strains = txt.count(taxonomy + "/")
        return num_conservation / num_all_strains

    def Funcs(self):
        return {a for b in [funcs(UniProtID) for UniProtID in self.UniProtID if funcs(UniProtID) is not None] for a in b}

    def shared_level(self):
        for rank in ranks:
            clades = taxon(rank)
            for clade in clades:
                if self.get_conservation(clade) > 0.98:
                    return rank
        return "Specific"


def list_core_gene(taxonomy, th_max, th_min):
    return {i for i in range(len_tsv) if th_min <= OrthoGroup(i).get_conservation(taxonomy) <= th_max}


def count_core_gene(taxonomy, th_max, th_min):
    count = 0
    list_cnsv = []
    for i in range(len_tsv):
        conservation = OrthoGroup(i).get_conservation(taxonomy)
        if th_min <= conservation <= th_max:
            count += 1
        list_cnsv.append(conservation)
    return count


def linage_specific_genes(taxonomy, th_max, th_min):
    rank = [rnk for rnk in ranks if taxonomy in get_taxonomy_list(rnk)][0]
    core_genes_true = set(list_core_gene(taxonomy, th_max, th_min))
    other_taxon = get_taxonomy_list(rank)
    other_taxon.remove(taxonomy)
    core_genes_others = [list_core_gene(other_taxa, th_max, th_min)
                         for other_taxa in other_taxon]
    core_genes_others = set([e for row in core_genes_others for e in row])
    linage_specific_genes = core_genes_true - core_genes_others
    return linage_specific_genes


def count_all_gene(species):
    for i in range(293):
        cmd = "cat ~/db/refseq/archaea/faa/seq{}.fasta | ".format(i) + \
              "head -n 1 | awk -F '@' '{print $2}'"
        name = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        name = name.stdout.decode("utf8").replace("\n", "")
        if name == species:
            cmd = "cat ~/db/refseq/archaea/faa/seq{}.fasta | ".format(i) + \
                  "grep '@' | wc -l"
            num = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
            num = num.stdout.decode("utf8")
            return int(num)


def numbers_to_func(numbers):
    funcs = [OrthoGroup(i).get_func() for i in numbers]
    funcs_c = [flat for flatten in [func[0]
                                    for func in funcs] for flat in flatten]
    funcs_f = [flat for flatten in [func[1]
                                    for func in funcs] for flat in flatten]
    funcs_p = [flat for flatten in [func[2]
                                    for func in funcs] for flat in flatten]
    df = pd.DataFrame(funcs_c)
    print(df.groupby(0).sum().sort_values(1, ascending=False))
    df = pd.DataFrame(funcs_f)
    print(df.groupby(0).sum().sort_values(1, ascending=False))
    df = pd.DataFrame(funcs_p)
    print(df.groupby(0).sum().sort_values(1, ascending=False))


def count_core_func(list_num):
    samples = [a for b in [OrthoGroup(ogn).Funcs()
                           for ogn in list_num] for a in b]
    res = collections.Counter(samples).most_common()
    return res


def Fisher_exact(clade_a, clade_b):
    Ar_set = list_core_gene(clade_a, 1, 0.98)
    Eu_set = list_core_gene(clade_b, 1, 0.98)
    lis_0 = count_core_func(Ar_set)
    lis_1 = count_core_func(Eu_set - Ar_set)
    df0 = pd.DataFrame(lis_0)
    df1 = pd.DataFrame(lis_1)
    for i in range(len(df0)):
        ar_gene = df0[1][i]
        try:
            oth_gene = int(df1[df1[0] == df0[0][i]][1])
        except TypeError:
            continue
        sum_ar = len(Ar_set)
        sum_oth = len(Eu_set) - len(Ar_set)
        ar_not_gene = sum_ar - ar_gene
        oth_not_gene = sum_oth - oth_gene
        matrix = [[ar_gene, ar_not_gene], [oth_gene, oth_not_gene]]
        oddsratio, pvalue = stats.fisher_exact(matrix)
        if pvalue < 0.1:
            sign = "ar>," if (
                (matrix[0][0]/matrix[0][1]) > (matrix[1][0]/matrix[1][1])) else "ar<,"
            print(sign, end="")
            print(df0[0][i].replace(",", ""), end=",")
            print(pvalue, end=",")
            print(str(matrix[0][0]) + "," + str(matrix[0][1]) +
                  "," + str(matrix[1][0]) + "," + str(matrix[1][1]))


def main1():
    for rank in ranks:
        taxonomy_list = [tax for tax in get_taxonomy_list(rank)]
        taxonomy_list = [tax for tax in taxonomy_list if tax in dedup]
        core_gene_num_list = [count_core_gene(
            taxonomy, 1, 0.98) for taxonomy in taxonomy_list]
        print(core_gene_num_list)


def main2():
    for i in range(23303):
        name = OrthoGroup(i).Name
        conservation = OrthoGroup(i).get_conservation("Archaea")
        if "ribosomal_protein" in name and conservation >= 0.98:
            print(name.split("_")[3])


def main3():
    numbers_to_func(list_core_gene("Euryarchaeota", 1, 0.97) -
                    list_core_gene("Archaea", 1, 0.97))


def main4():
    lsg_list = linage_specific_genes("Pyrococcus", 1, 0.98)
    print([OrthoGroup(i).Name for i in lsg_list])
    print(len(lsg_list))


def main5():
    li = list(set(list_core_gene("Crenarchaeota", 1, 0.98)) -
              set(list_core_gene("Archaea", 1, 0.98)))
    li = [OrthoGroup(i).get_func()
          for i in li if OrthoGroup(i).get_func() is not None]
    li = [a for b in li for a in b]
    di = collections.Counter(li)
    df = pd.DataFrame(di.most_common())
    df.to_csv("Cren_core.csv", header=False, index=False)


def decrease():
    res = [str(len(list_core_gene(taxa, 1, 0.98)))
           for taxa in taxon("Species")]
    print("\n".join(res))


def prevarence():
    ranks.append("Specific")
    for rank in ranks:
        print(rank)
        res = [[OrthoGroup(i).get_conservation("Archaea"), OrthoGroup(
            i).shared_level(), 1] for i in range(len_tsv)]
        df = pd.DataFrame(res)
        df = df[df[1] == rank].drop(1, axis=1)
        df1 = df.groupby([0]).count()
        print("\n".join([",".join([str(i/292), str(df1[df1.index == i/292].values[0][0])]) if i /
                         292 in df1.index else ",".join([str(i/292), str(0)]) for i in range(1, 293)]))


if __name__ == "__main__":
#    core = [0, 1, 2, 4, 5, 6, 518, 521, 522, 12, 525, 15, 527, 18, 531, 532, 21, 22, 535, 25, 537, 539, 542, 33, 547, 37, 551, 41, 42, 48, 49, 50, 51, 53, 60, 61, 69, 71, 78, 85, 87, 88, 90, 91, 95, 103, 106, 107, 109, 110, 112, 116, 118, 123, 124, 128, 129, 130, 132, 139, 140, 143, 148, 151, 153, 156, 158, 159, 169, 171, 176, 188, 192, 195, 200, 202, 206, 209, 211, 217, 218, 219, 221, 228, 230, 231, 234, 236, 240, 241, 245, 247, 248, 249, 252, 253, 257, 259, 260, 261, 264, 267, 268, 269, 270, 272, 273, 274, 275, 276, 281, 282,
#            284, 285, 286, 288, 289, 290, 291, 292, 294, 296, 298, 299, 300, 301, 302, 304, 305, 306, 307, 308, 309, 311, 313, 316, 318, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 334, 335, 338, 339, 342, 343, 344, 347, 348, 349, 350, 353, 355, 356, 357, 358, 359, 362, 365, 372, 373, 375, 378, 379, 381, 386, 388, 390, 391, 394, 396, 397, 398, 401, 408, 411, 412, 414, 415, 419, 420, 423, 424, 427, 430, 432, 435, 437, 441, 445, 448, 450, 452, 460, 461, 462, 463, 466, 467, 472, 480, 481, 483, 487, 488, 493, 497, 499, 506, 511]
    for i in range(len_tsv):
        if "ribosomal_protein" in OrthoGroup(i).Name and OrthoGroup(i).get_conservation("Archaea") > 0.98:
            print(OrthoGroup(i).Name)
