import re


def replace_string(string):
    string = re.sub("\"", "", string)
    string = re.sub("\t", "", string)
    string = re.sub("\n", "", string)
    string = re.sub("'", "", string)
    string = re.sub("\(", "_", string)
    string = re.sub("\)", "_", string)
    string = re.sub("\|", "_", string)
    string = re.sub("\[", "", string)
    string = re.sub(",", "", string)
    string = re.sub("\.", "", string)
    string = re.sub("@", "", string)
    string = re.sub("\*", "", string)
    string = re.sub("]", "", string)
    string = re.sub("\+", "", string)
    string = re.sub("{", "", string)
    string = re.sub("}", "", string)
    string = re.sub("\^", "", string)
    string = re.sub("/", "", string)
    string = re.sub("\$", "", string)
    string = re.sub(" ", "_", string)
    string = re.sub("\\\\", "", string)
    string = re.sub(";", "", string)
    string = re.sub("/", "", string)
    string = re.sub("<", "", string)
    string = re.sub(":", "", string)
    string = re.sub("¥", "", string)
    string = re.sub(" ", "_", string)
    return string


f = open("swiss_prot_archaea_2.fasta", "r")
lines = f.readlines()
f.close()

text = ""
for line in lines:
    l = replace_string(line)
    text = text + l + "\n"

print(text)
