import pandas as pd


df = pd.read_csv("UniProt_GO.tsv", sep="\t", header=None)


def funcs(UniProtID):
    try:
        return str(df[df[0] == UniProtID][1].values[0]).split("@")
    except IndexError:
        return None
