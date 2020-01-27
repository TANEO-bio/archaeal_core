from ast import literal_eval


f = open("Ortho_Func.tsv")
lines = f.readlines()
f.close()


def process(func_foo):
    func_foo = func_foo.replace("\n", "").replace(
        "C:", "").replace("F:", "").replace("P:", "")
    if func_foo != "":
        func_foo = literal_eval(func_foo)
        if type(func_foo[1]) is tuple:
            func_foo_sum = sum([i[1] for i in func_foo])
            func_foo = [(func[0], func[1]/func_foo_sum) for func in func_foo]
        else:
            func_foo = [(func_foo[0], 1)]
    else:
        func_foo = []
    return func_foo


def num_to_func(num):
    funcs = lines[num].split("\t")
    funcs = [funcs[1], funcs[2], funcs[3]]
    return list(map(process, funcs))
