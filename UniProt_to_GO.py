f = open("/home/t18476nt/db/GO/test.txt")
lines = f.readlines()
f.close()

len_f = len(lines)


class term:
    def __init__(self, term):
        self.term = term
        self.pos = [i for i in range(
            len_f) if "id: {}".format(term) in lines[i]][0]
        self.name = lines[self.pos].split("name: ")[1].split("\t")[0]
        tmp = lines[self.pos].split("\t")
        if len(tmp) >= 3:
            self.anc = {"GO:" + elm.split("GO:")[1].split(" !")[0]
                        for elm in tmp[2:]}
            if self.term in self.anc:
                self.anc.remove(self.term)
        else:
            self.anc = None

    def tree(self):
        funcs = [self.name]
        query = set(self.anc)
        searched = set()
        while len(query) > 0:
            ins = term(query.pop())
            funcs.append(ins.name.replace("\n", ""))
            if ins.anc is not None:
                if ins.term not in searched:
                    query = query | ins.anc
                    searched.add(ins.term)
        return funcs


