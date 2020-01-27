from Taxid_Taxonomy import get_taxonomy_list


h = open("dedup.txt")
dedup = [string.replace("\n", "") for string in h.readlines()]
h.close()

ranks = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

dic = {el[0]: el[1] for el in [ele for ele in [(taxa, rank) * (taxa in get_taxonomy_list(rank))
                                               for rank in ranks for taxa in dedup] if ele != ()]}


def taxon(rank):
    return [k for k, v in dic.items() if v == rank]
