#!usr/bin/env python
import glob
import sys
import subprocess
from ete3 import NCBITaxa
from Bio import SeqIO
import re


def _taxid2lineage(taxid):
    lineage = ncbi.get_lineage(taxid)
    lineage_dict = dict()
    names = ncbi.get_taxid_translator(lineage)
    for rank in ranks:
        for k, v in ncbi.get_rank(lineage).items():
            if v == rank:
                lineage_dict.update({v: names[k]})
    return lineage_dict


ranks = ['superkingdom', 'phylum', 'class',
         'order', 'family', 'genus', 'species']
ncbi = NCBITaxa()
arg = int(sys.argv[1])
path_list = glob.glob("/home/t18476nt/db/refseq/archaea/gbff2/*.gbff")
path_list = path_list[arg::32]
for path in path_list:
    print(path)
    key = 0
    text = ""
    for rec in SeqIO.parse(path, "genbank"):
        accession = rec.id
        source = rec.annotations['source']
        for feat in rec.features:
            if "pseudo" not in feat.qualifiers and "pseudogene" not in feat.qualifiers:
                if feat.qualifiers.get("db_xref") is not None:
                    dbxref = feat.qualifiers.get("db_xref")
                    tax_id = [tax for tax in dbxref if 'taxon:' in tax]
                    if tax_id != []:
                        taxid = re.sub(r'\D', '', tax_id[0])
                if feat.qualifiers.get("db_xref") is not None:
                    dbxref = feat.qualifiers.get("db_xref")
                    geneid = [gene for gene in dbxref if 'GeneID:' in gene]
                    if geneid != []:
                        gene_id = re.sub(r'\D', '', geneid[0])
                if feat.qualifiers.get("protein_id") is not None:
                    protein_id = feat.qualifiers.get("protein_id")[0]
                if feat.qualifiers.get("product") is not None:
                    product = feat.qualifiers.get("product")[0]
                if feat.location.start is not None:
                    start = str(feat.location.start)
                if feat.location.end is not None:
                    end = str(feat.location.end)
                if feat.qualifiers.get("translation") is not None:
                    seq = feat.qualifiers.get("translation")[0]
                    ids = [accession, source, taxid, gene_id,
                           protein_id, product, start, end]
                    new_ids = []
                    for j in ids:
                        j = re.sub("\"", "_", j)
                        j = re.sub("'", "_", j)
                        j = re.sub("\(", "_", j)
                        j = re.sub("\)", "_", j)
                        j = re.sub("\|", "_", j)
                        j = re.sub("\[", "_", j)
                        j = re.sub(",", "_", j)
                        j = re.sub("\.", "_", j)
                        j = re.sub("@", "_", j)
                        j = re.sub("\*", "_", j)
                        j = re.sub("]", "_", j)
                        j = re.sub("\+", "_", j)
                        j = re.sub("{", "_", j)
                        j = re.sub("}", "_", j)
                        j = re.sub("\^", "_", j)
                        j = re.sub("/", "_", j)
                        j = re.sub("\$", "_", j)
                        j = re.sub(" ", "_", j)
                        j = re.sub("\\\\", "_", j)
                        j = re.sub(";", "_", j)
                        j = re.sub("/", "_", j)
                        j = re.sub("<", "_", j)
                        j = re.sub(">", "_", j)
                        j = re.sub(":", "_", j)
                        j = re.sub("¥", "_", j)
                        new_ids.append(j)
                    row = ">{0}@{1}@{2}@{3}@{4}@{5}@{6}_{7}".format(
                        new_ids[0], new_ids[1], new_ids[2], new_ids[3], new_ids[4],
                        new_ids[5], new_ids[6], new_ids[7]) + "\n" + seq + "\n"
                    text = text + row
                    try:
                        taxonomy = _taxid2lineage(taxid)
                        key = 1

                    except:
                        subprocess.run("rm " + path, shell=True)

        if key == 1:
            path_w = '/home/t18476nt/db/refseq/archaea/faa2/seq{}.fasta'.format(arg)
            f = open(path_w, mode='w')
            f.write(text)
            f.close
            print(arg)
            arg = arg + 32
            break
