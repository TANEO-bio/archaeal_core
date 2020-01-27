import urllib.request
import sys
from bs4 import BeautifulSoup
from UniProt_to_GO import term
import subprocess


def uniprotid_to_go(uniprotid):
    if uniprotid is not None:
        url = "https://www.uniprot.org/uniprot/{}.xml".format(uniprotid)
        xml = urllib.request.urlopen(url, timeout=10)
        soup = BeautifulSoup(xml, 'lxml')
        tags = soup.find_all("property", type="term")
        GO = []
        for tag in tags:
            func = tag["value"]
            GO.append(func)
        return GO
    else:
        return None


def uniprotid_to_GO_Term(uniprotid):
    if uniprotid is not None:
        url = "https://www.uniprot.org/uniprot/{}.xml".format(uniprotid)
        xml = urllib.request.urlopen(url, timeout=10)
        soup = BeautifulSoup(xml, 'lxml')
        tags = soup.find_all("dbreference", type="GO")
        GO = []
        for tag in tags:
            func = tag["id"]
            GO.append(func)
        return GO
    else:
        return None

num = int(sys.argv[1])

f = open("/home/t18476nt/db/swiss_prot_archaea/swiss_prot_archaea_2.fasta")
lines = f.readlines()
f.close()

txt = ""
p = num
for line in lines[num::50]:
    if ">" in line:
        id = line.split("_")[1]
        try:
            GOs = uniprotid_to_GO_Term(id)
        except:
            GOs = []
        terms = [term(GO).tree() for GO in GOs]
        terms = {a for b in terms for a in b}
        lin = id + "\t" + "@".join(terms) + "\n"
        txt = txt + lin
        subprocess.run("echo {} > progress{}.txt".format(p, num), shell=True)
        p = p + 15

g = open("UniProt_GO{}.tsv".format(num), "w")
g.write(txt)
g.close()
