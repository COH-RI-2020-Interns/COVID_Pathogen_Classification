import itertools

d = {}
l = []
contd = 0;


def retrieve_RNA(seq):
    for AA in seq:
        with open("ct.txt") as fin:
            for line in fin:
                items = line.split()
                p = items[2:]
                if str(AA) == str(items[0]):
                    l.append(p)

            posbl = list(itertools.product(*l))

            print ("****** All Possible Sequences ******\n")
            for sequence in posbl:
                for codon in sequence:
                    print(str(codon) + " "),
                print ("\n ")
            return posbl


retrieve_RNA('LPEETMDFMVQHTSFKEMKKNPMTNYTTVPQELMDHSISPFMRKGM')


x = [('UUA', 'CUA', 'UUA', 'CUA', 'UUA', 'CUA'),
 ('UUA', 'CUA', 'UUA', 'CUA', 'UUA', 'CUU'),
 ('UUA', 'CUA', 'UUA', 'CUA', 'UUA', 'CUG'),
 ('UUA', 'CUA', 'UUA', 'CUA', 'UUA', 'CUC'),
 ('UUA', 'CUA', 'UUA', 'CUA', 'UUG', 'CUA'),
 ('UUA', 'CUA', 'UUA', 'CUA', 'UUG', 'CUU'),
 ('UUA', 'CUA', 'UUA', 'CUA', 'UUG', 'CUG'),
 ('UUA', 'CUA', 'UUA', 'CUA', 'UUG', 'CUC'),
 ('UUA', 'CUA', 'UUA', 'CUU', 'UUA', 'CUA'),
 ('UUA', 'CUA', 'UUA', 'CUU', 'UUA', 'CUU'),
 ('UUA', 'CUA', 'UUA', 'CUU', 'UUA', 'CUG'),
 ('UUA', 'CUA', 'UUA', 'CUU', 'UUA', 'CUC'),
 ('UUA', 'CUA', 'UUA', 'CUU', 'UUG', 'CUA'),
 ('UUA', 'CUA', 'UUA', 'CUU', 'UUG', 'CUU'),
 ('UUA', 'CUA', 'UUA', 'CUU', 'UUG', 'CUG'),
 ('UUA', 'CUA', 'UUA', 'CUU', 'UUG', 'CUC'),
 ('UUA', 'CUA', 'UUA', 'CUG', 'UUA', 'CUA'),
 ('UUA', 'CUA', 'UUA', 'CUG', 'UUA', 'CUU'),
 ('UUA', 'CUA', 'UUA', 'CUG', 'UUA', 'CUG'),
 ('UUA', 'CUA', 'UUA', 'CUG', 'UUA', 'CUC'),
 ('UUA', 'CUA', 'UUA', 'CUG', 'UUG', 'CUA'),
 ('UUA', 'CUA', 'UUA', 'CUG', 'UUG', 'CUU'),
 ('UUA', 'CUA', 'UUA', 'CUG', 'UUG', 'CUG'),
 ('UUA', 'CUA', 'UUA', 'CUG', 'UUG', 'CUC'),
 ('UUA', 'CUA', 'UUA', 'CUC', 'UUA', 'CUA'),
 ('UUA', 'CUA', 'UUA', 'CUC', 'UUA', 'CUU'),
 ('UUA', 'CUA', 'UUA', 'CUC', 'UUA', 'CUG'),
 ('UUA', 'CUA', 'UUA', 'CUC', 'UUA', 'CUC'),
 ('UUA', 'CUA', 'UUA', 'CUC', 'UUG', 'CUA'),
 ('UUA', 'CUA', 'UUA', 'CUC', 'UUG', 'CUU'),
 ('UUA', 'CUA', 'UUA', 'CUC', 'UUG', 'CUG'),
 ('UUA', 'CUA', 'UUA', 'CUC', 'UUG', 'CUC'),
 ('UUA', 'CUA', 'UUG', 'CUA', 'UUA', 'CUA'),
 ('UUA', 'CUA', 'UUG', 'CUA', 'UUA', 'CUU'),
 ('UUA', 'CUA', 'UUG', 'CUA', 'UUA', 'CUG'),
 ('UUA', 'CUA', 'UUG', 'CUA', 'UUA', 'CUC'),
 ('UUA', 'CUA', 'UUG', 'CUA', 'UUG', 'CUA'),
 ('UUA', 'CUA', 'UUG', 'CUA', 'UUG', 'CUU'),
 ('UUA', 'CUA', 'UUG', 'CUA', 'UUG', 'CUG'),
 ('UUA', 'CUA', 'UUG', 'CUA', 'UUG', 'CUC'),
 ('UUA', 'CUA', 'UUG', 'CUU', 'UUA', 'CUA'),
 ('UUA', 'CUA', 'UUG', 'CUU', 'UUA', 'CUU'),
 ('UUA', 'CUA', 'UUG', 'CUU', 'UUA', 'CUG'),
 ('UUA', 'CUA', 'UUG', 'CUU', 'UUA', 'CUC'),
 ('UUA', 'CUA', 'UUG', 'CUU', 'UUG', 'CUA'),
 ('UUA', 'CUA', 'UUG', 'CUU', 'UUG', 'CUU'),
 ('UUA', 'CUA', 'UUG', 'CUU', 'UUG', 'CUG'),
 ('UUA', 'CUA', 'UUG', 'CUU', 'UUG', 'CUC'),
 ('UUA', 'CUA', 'UUG', 'CUG', 'UUA', 'CUA'),
 ('UUA', 'CUA', 'UUG', 'CUG', 'UUA', 'CUU'),
 ('UUA', 'CUA', 'UUG', 'CUG', 'UUA', 'CUG'),
 ('UUA', 'CUA', 'UUG', 'CUG', 'UUA', 'CUC'),
 ('UUA', 'CUA', 'UUG', 'CUG', 'UUG', 'CUA'),
 ('UUA', 'CUA', 'UUG', 'CUG', 'UUG', 'CUU'),
 ('UUA', 'CUA', 'UUG', 'CUG', 'UUG', 'CUG'),
 ('UUA', 'CUA', 'UUG', 'CUG', 'UUG', 'CUC'),
 ('UUA', 'CUA', 'UUG', 'CUC', 'UUA', 'CUA'),
 ('UUA', 'CUA', 'UUG', 'CUC', 'UUA', 'CUU'),
 ('UUA', 'CUA', 'UUG', 'CUC', 'UUA', 'CUG'),
 ('UUA', 'CUA', 'UUG', 'CUC', 'UUA', 'CUC'),
 ('UUA', 'CUA', 'UUG', 'CUC', 'UUG', 'CUA'),
 ('UUA', 'CUA', 'UUG', 'CUC', 'UUG', 'CUU'),
 ('UUA', 'CUA', 'UUG', 'CUC', 'UUG', 'CUG'),
 ('UUA', 'CUA', 'UUG', 'CUC', 'UUG', 'CUC'),
 ('UUA', 'CUU', 'UUA', 'CUA', 'UUA', 'CUA'),
 ('UUA', 'CUU', 'UUA', 'CUA', 'UUA', 'CUU'),
 ('UUA', 'CUU', 'UUA', 'CUA', 'UUA', 'CUG'),
 ('UUA', 'CUU', 'UUA', 'CUA', 'UUA', 'CUC'),
 ('UUA', 'CUU', 'UUA', 'CUA', 'UUG', 'CUA'),
 ('UUA', 'CUU', 'UUA', 'CUA', 'UUG', 'CUU'),
 ('UUA', 'CUU', 'UUA', 'CUA', 'UUG', 'CUG'),
 ('UUA', 'CUU', 'UUA', 'CUA', 'UUG', 'CUC'),
 ('UUA', 'CUU', 'UUA', 'CUU', 'UUA', 'CUA'),
 ('UUA', 'CUU', 'UUA', 'CUU', 'UUA', 'CUU'),
 ('UUA', 'CUU', 'UUA', 'CUU', 'UUA', 'CUG'),
 ('UUA', 'CUU', 'UUA', 'CUU', 'UUA', 'CUC'),
 ('UUA', 'CUU', 'UUA', 'CUU', 'UUG', 'CUA'),
 ('UUA', 'CUU', 'UUA', 'CUU', 'UUG', 'CUU'),
 ('UUA', 'CUU', 'UUA', 'CUU', 'UUG', 'CUG'),
 ('UUA', 'CUU', 'UUA', 'CUU', 'UUG', 'CUC'),
 ('UUA', 'CUU', 'UUA', 'CUG', 'UUA', 'CUA'),
 ('UUA', 'CUU', 'UUA', 'CUG', 'UUA', 'CUU'),
 ('UUA', 'CUU', 'UUA', 'CUG', 'UUA', 'CUG'),
 ('UUA', 'CUU', 'UUA', 'CUG', 'UUA', 'CUC'),
 ('UUA', 'CUU', 'UUA', 'CUG', 'UUG', 'CUA'),
 ('UUA', 'CUU', 'UUA', 'CUG', 'UUG', 'CUU'),
 ('UUA', 'CUU', 'UUA', 'CUG', 'UUG', 'CUG'),
 ('UUA', 'CUU', 'UUA', 'CUG', 'UUG', 'CUC'),
 ('UUA', 'CUU', 'UUA', 'CUC', 'UUA', 'CUA'),
 ('UUA', 'CUU', 'UUA', 'CUC', 'UUA', 'CUU'),
 ('UUA', 'CUU', 'UUA', 'CUC', 'UUA', 'CUG'),
 ('UUA', 'CUU', 'UUA', 'CUC', 'UUA', 'CUC'),
 ('UUA', 'CUU', 'UUA', 'CUC', 'UUG', 'CUA'),
 ('UUA', 'CUU', 'UUA', 'CUC', 'UUG', 'CUU'),
 ('UUA', 'CUU', 'UUA', 'CUC', 'UUG', 'CUG'),
 ('UUA', 'CUU', 'UUA', 'CUC', 'UUG', 'CUC'),
 ('UUA', 'CUU', 'UUG', 'CUA', 'UUA', 'CUA'),
 ('UUA', 'CUU', 'UUG', 'CUA', 'UUA', 'CUU'),
 ('UUA', 'CUU', 'UUG', 'CUA', 'UUA', 'CUG'),
 ('UUA', 'CUU', 'UUG', 'CUA', 'UUA', 'CUC'),
 ('UUA', 'CUU', 'UUG', 'CUA', 'UUG', 'CUA'),
 ('UUA', 'CUU', 'UUG', 'CUA', 'UUG', 'CUU'),
 ('UUA', 'CUU', 'UUG', 'CUA', 'UUG', 'CUG'),
 ('UUA', 'CUU', 'UUG', 'CUA', 'UUG', 'CUC'),
 ('UUA', 'CUU', 'UUG', 'CUU', 'UUA', 'CUA'),
 ('UUA', 'CUU', 'UUG', 'CUU', 'UUA', 'CUU'),
 ('UUA', 'CUU', 'UUG', 'CUU', 'UUA', 'CUG'),
 ('UUA', 'CUU', 'UUG', 'CUU', 'UUA', 'CUC'),
 ('UUA', 'CUU', 'UUG', 'CUU', 'UUG', 'CUA'),
 ('UUA', 'CUU', 'UUG', 'CUU', 'UUG', 'CUU'),
 ('UUA', 'CUU', 'UUG', 'CUU', 'UUG', 'CUG'),
 ('UUA', 'CUU', 'UUG', 'CUU', 'UUG', 'CUC'),
 ('UUA', 'CUU', 'UUG', 'CUG', 'UUA', 'CUA'),
 ('UUA', 'CUU', 'UUG', 'CUG', 'UUA', 'CUU'),
 ('UUA', 'CUU', 'UUG', 'CUG', 'UUA', 'CUG'),
 ('UUA', 'CUU', 'UUG', 'CUG', 'UUA', 'CUC'),
 ('UUA', 'CUU', 'UUG', 'CUG', 'UUG', 'CUA'),
 ('UUA', 'CUU', 'UUG', 'CUG', 'UUG', 'CUU'),
 ('UUA', 'CUU', 'UUG', 'CUG', 'UUG', 'CUG'),
 ('UUA', 'CUU', 'UUG', 'CUG', 'UUG', 'CUC'),
 ('UUA', 'CUU', 'UUG', 'CUC', 'UUA', 'CUA'),
 ('UUA', 'CUU', 'UUG', 'CUC', 'UUA', 'CUU'),
 ('UUA', 'CUU', 'UUG', 'CUC', 'UUA', 'CUG'),
 ('UUA', 'CUU', 'UUG', 'CUC', 'UUA', 'CUC'),
 ('UUA', 'CUU', 'UUG', 'CUC', 'UUG', 'CUA'),
 ('UUA', 'CUU', 'UUG', 'CUC', 'UUG', 'CUU'),
 ('UUA', 'CUU', 'UUG', 'CUC', 'UUG', 'CUG'),
 ('UUA', 'CUU', 'UUG', 'CUC', 'UUG', 'CUC'),
 ('UUA', 'CUG', 'UUA', 'CUA', 'UUA', 'CUA'),
 ('UUA', 'CUG', 'UUA', 'CUA', 'UUA', 'CUU'),
 ('UUA', 'CUG', 'UUA', 'CUA', 'UUA', 'CUG'),
 ('UUA', 'CUG', 'UUA', 'CUA', 'UUA', 'CUC'),
 ('UUA', 'CUG', 'UUA', 'CUA', 'UUG', 'CUA'),
 ('UUA', 'CUG', 'UUA', 'CUA', 'UUG', 'CUU'),
 ('UUA', 'CUG', 'UUA', 'CUA', 'UUG', 'CUG'),
 ('UUA', 'CUG', 'UUA', 'CUA', 'UUG', 'CUC'),
 ('UUA', 'CUG', 'UUA', 'CUU', 'UUA', 'CUA'),
 ('UUA', 'CUG', 'UUA', 'CUU', 'UUA', 'CUU'),
 ('UUA', 'CUG', 'UUA', 'CUU', 'UUA', 'CUG'),
 ('UUA', 'CUG', 'UUA', 'CUU', 'UUA', 'CUC'),
 ('UUA', 'CUG', 'UUA', 'CUU', 'UUG', 'CUA'),
 ('UUA', 'CUG', 'UUA', 'CUU', 'UUG', 'CUU'),
 ('UUA', 'CUG', 'UUA', 'CUU', 'UUG', 'CUG'),
 ('UUA', 'CUG', 'UUA', 'CUU', 'UUG', 'CUC'),
 ('UUA', 'CUG', 'UUA', 'CUG', 'UUA', 'CUA'),
 ('UUA', 'CUG', 'UUA', 'CUG', 'UUA', 'CUU'),
 ('UUA', 'CUG', 'UUA', 'CUG', 'UUA', 'CUG'),
 ('UUA', 'CUG', 'UUA', 'CUG', 'UUA', 'CUC'),
 ('UUA', 'CUG', 'UUA', 'CUG', 'UUG', 'CUA'),
 ('UUA', 'CUG', 'UUA', 'CUG', 'UUG', 'CUU'),
 ('UUA', 'CUG', 'UUA', 'CUG', 'UUG', 'CUG'),
 ('UUA', 'CUG', 'UUA', 'CUG', 'UUG', 'CUC'),
 ('UUA', 'CUG', 'UUA', 'CUC', 'UUA', 'CUA'),
 ('UUA', 'CUG', 'UUA', 'CUC', 'UUA', 'CUU'),
 ('UUA', 'CUG', 'UUA', 'CUC', 'UUA', 'CUG'),
 ('UUA', 'CUG', 'UUA', 'CUC', 'UUA', 'CUC'),
 ('UUA', 'CUG', 'UUA', 'CUC', 'UUG', 'CUA'),
 ('UUA', 'CUG', 'UUA', 'CUC', 'UUG', 'CUU'),
 ('UUA', 'CUG', 'UUA', 'CUC', 'UUG', 'CUG'),
 ('UUA', 'CUG', 'UUA', 'CUC', 'UUG', 'CUC'),
 ('UUA', 'CUG', 'UUG', 'CUA', 'UUA', 'CUA'),
 ('UUA', 'CUG', 'UUG', 'CUA', 'UUA', 'CUU'),
 ('UUA', 'CUG', 'UUG', 'CUA', 'UUA', 'CUG'),
 ('UUA', 'CUG', 'UUG', 'CUA', 'UUA', 'CUC'),
 ('UUA', 'CUG', 'UUG', 'CUA', 'UUG', 'CUA'),
 ('UUA', 'CUG', 'UUG', 'CUA', 'UUG', 'CUU'),
 ('UUA', 'CUG', 'UUG', 'CUA', 'UUG', 'CUG'),
 ('UUA', 'CUG', 'UUG', 'CUA', 'UUG', 'CUC'),
 ('UUA', 'CUG', 'UUG', 'CUU', 'UUA', 'CUA'),
 ('UUA', 'CUG', 'UUG', 'CUU', 'UUA', 'CUU'),
 ('UUA', 'CUG', 'UUG', 'CUU', 'UUA', 'CUG'),
 ('UUA', 'CUG', 'UUG', 'CUU', 'UUA', 'CUC'),
 ('UUA', 'CUG', 'UUG', 'CUU', 'UUG', 'CUA'),
 ('UUA', 'CUG', 'UUG', 'CUU', 'UUG', 'CUU'),
 ('UUA', 'CUG', 'UUG', 'CUU', 'UUG', 'CUG'),
 ('UUA', 'CUG', 'UUG', 'CUU', 'UUG', 'CUC'),
 ('UUA', 'CUG', 'UUG', 'CUG', 'UUA', 'CUA'),
 ('UUA', 'CUG', 'UUG', 'CUG', 'UUA', 'CUU'),
 ('UUA', 'CUG', 'UUG', 'CUG', 'UUA', 'CUG'),
 ('UUA', 'CUG', 'UUG', 'CUG', 'UUA', 'CUC'),
 ('UUA', 'CUG', 'UUG', 'CUG', 'UUG', 'CUA'),
 ('UUA', 'CUG', 'UUG', 'CUG', 'UUG', 'CUU'),
 ('UUA', 'CUG', 'UUG', 'CUG', 'UUG', 'CUG'),
 ('UUA', 'CUG', 'UUG', 'CUG', 'UUG', 'CUC'),
 ('UUA', 'CUG', 'UUG', 'CUC', 'UUA', 'CUA'),
 ('UUA', 'CUG', 'UUG', 'CUC', 'UUA', 'CUU'),
 ('UUA', 'CUG', 'UUG', 'CUC', 'UUA', 'CUG'),
 ('UUA', 'CUG', 'UUG', 'CUC', 'UUA', 'CUC'),
 ('UUA', 'CUG', 'UUG', 'CUC', 'UUG', 'CUA'),
 ('UUA', 'CUG', 'UUG', 'CUC', 'UUG', 'CUU'),
 ('UUA', 'CUG', 'UUG', 'CUC', 'UUG', 'CUG'),
 ('UUA', 'CUG', 'UUG', 'CUC', 'UUG', 'CUC'),
 ('UUA', 'CUC', 'UUA', 'CUA', 'UUA', 'CUA'),
 ('UUA', 'CUC', 'UUA', 'CUA', 'UUA', 'CUU'),
 ('UUA', 'CUC', 'UUA', 'CUA', 'UUA', 'CUG'),
 ('UUA', 'CUC', 'UUA', 'CUA', 'UUA', 'CUC'),
 ('UUA', 'CUC', 'UUA', 'CUA', 'UUG', 'CUA'),
 ('UUA', 'CUC', 'UUA', 'CUA', 'UUG', 'CUU'),
 ('UUA', 'CUC', 'UUA', 'CUA', 'UUG', 'CUG'),
 ('UUA', 'CUC', 'UUA', 'CUA', 'UUG', 'CUC'),
 ('UUA', 'CUC', 'UUA', 'CUU', 'UUA', 'CUA'),
 ('UUA', 'CUC', 'UUA', 'CUU', 'UUA', 'CUU'),
 ('UUA', 'CUC', 'UUA', 'CUU', 'UUA', 'CUG'),
 ('UUA', 'CUC', 'UUA', 'CUU', 'UUA', 'CUC'),
 ('UUA', 'CUC', 'UUA', 'CUU', 'UUG', 'CUA'),
 ('UUA', 'CUC', 'UUA', 'CUU', 'UUG', 'CUU'),
 ('UUA', 'CUC', 'UUA', 'CUU', 'UUG', 'CUG'),
 ('UUA', 'CUC', 'UUA', 'CUU', 'UUG', 'CUC'),
 ('UUA', 'CUC', 'UUA', 'CUG', 'UUA', 'CUA'),
 ('UUA', 'CUC', 'UUA', 'CUG', 'UUA', 'CUU'),
 ('UUA', 'CUC', 'UUA', 'CUG', 'UUA', 'CUG'),
 ('UUA', 'CUC', 'UUA', 'CUG', 'UUA', 'CUC'),
 ('UUA', 'CUC', 'UUA', 'CUG', 'UUG', 'CUA'),
 ('UUA', 'CUC', 'UUA', 'CUG', 'UUG', 'CUU'),
 ('UUA', 'CUC', 'UUA', 'CUG', 'UUG', 'CUG'),
 ('UUA', 'CUC', 'UUA', 'CUG', 'UUG', 'CUC'),
 ('UUA', 'CUC', 'UUA', 'CUC', 'UUA', 'CUA'),
 ('UUA', 'CUC', 'UUA', 'CUC', 'UUA', 'CUU'),
 ('UUA', 'CUC', 'UUA', 'CUC', 'UUA', 'CUG'),
 ('UUA', 'CUC', 'UUA', 'CUC', 'UUA', 'CUC'),
 ('UUA', 'CUC', 'UUA', 'CUC', 'UUG', 'CUA'),
 ('UUA', 'CUC', 'UUA', 'CUC', 'UUG', 'CUU'),
 ('UUA', 'CUC', 'UUA', 'CUC', 'UUG', 'CUG'),
 ('UUA', 'CUC', 'UUA', 'CUC', 'UUG', 'CUC'),
 ('UUA', 'CUC', 'UUG', 'CUA', 'UUA', 'CUA'),
 ('UUA', 'CUC', 'UUG', 'CUA', 'UUA', 'CUU'),
 ('UUA', 'CUC', 'UUG', 'CUA', 'UUA', 'CUG'),
 ('UUA', 'CUC', 'UUG', 'CUA', 'UUA', 'CUC'),
 ('UUA', 'CUC', 'UUG', 'CUA', 'UUG', 'CUA'),
 ('UUA', 'CUC', 'UUG', 'CUA', 'UUG', 'CUU'),
 ('UUA', 'CUC', 'UUG', 'CUA', 'UUG', 'CUG'),
 ('UUA', 'CUC', 'UUG', 'CUA', 'UUG', 'CUC'),
 ('UUA', 'CUC', 'UUG', 'CUU', 'UUA', 'CUA'),
 ('UUA', 'CUC', 'UUG', 'CUU', 'UUA', 'CUU'),
 ('UUA', 'CUC', 'UUG', 'CUU', 'UUA', 'CUG'),
 ('UUA', 'CUC', 'UUG', 'CUU', 'UUA', 'CUC'),
 ('UUA', 'CUC', 'UUG', 'CUU', 'UUG', 'CUA'),
 ('UUA', 'CUC', 'UUG', 'CUU', 'UUG', 'CUU'),
 ('UUA', 'CUC', 'UUG', 'CUU', 'UUG', 'CUG'),
 ('UUA', 'CUC', 'UUG', 'CUU', 'UUG', 'CUC'),
 ('UUA', 'CUC', 'UUG', 'CUG', 'UUA', 'CUA'),
 ('UUA', 'CUC', 'UUG', 'CUG', 'UUA', 'CUU'),
 ('UUA', 'CUC', 'UUG', 'CUG', 'UUA', 'CUG'),
 ('UUA', 'CUC', 'UUG', 'CUG', 'UUA', 'CUC'),
 ('UUA', 'CUC', 'UUG', 'CUG', 'UUG', 'CUA'),
 ('UUA', 'CUC', 'UUG', 'CUG', 'UUG', 'CUU'),
 ('UUA', 'CUC', 'UUG', 'CUG', 'UUG', 'CUG'),
 ('UUA', 'CUC', 'UUG', 'CUG', 'UUG', 'CUC'),
 ('UUA', 'CUC', 'UUG', 'CUC', 'UUA', 'CUA'),
 ('UUA', 'CUC', 'UUG', 'CUC', 'UUA', 'CUU'),
 ('UUA', 'CUC', 'UUG', 'CUC', 'UUA', 'CUG'),
 ('UUA', 'CUC', 'UUG', 'CUC', 'UUA', 'CUC'),
 ('UUA', 'CUC', 'UUG', 'CUC', 'UUG', 'CUA'),
 ('UUA', 'CUC', 'UUG', 'CUC', 'UUG', 'CUU'),
 ('UUA', 'CUC', 'UUG', 'CUC', 'UUG', 'CUG'),
 ('UUA', 'CUC', 'UUG', 'CUC', 'UUG', 'CUC'),
 ('UUG', 'CUA', 'UUA', 'CUA', 'UUA', 'CUA'),
 ('UUG', 'CUA', 'UUA', 'CUA', 'UUA', 'CUU'),
 ('UUG', 'CUA', 'UUA', 'CUA', 'UUA', 'CUG'),
 ('UUG', 'CUA', 'UUA', 'CUA', 'UUA', 'CUC'),
 ('UUG', 'CUA', 'UUA', 'CUA', 'UUG', 'CUA'),
 ('UUG', 'CUA', 'UUA', 'CUA', 'UUG', 'CUU'),
 ('UUG', 'CUA', 'UUA', 'CUA', 'UUG', 'CUG'),
 ('UUG', 'CUA', 'UUA', 'CUA', 'UUG', 'CUC'),
 ('UUG', 'CUA', 'UUA', 'CUU', 'UUA', 'CUA'),
 ('UUG', 'CUA', 'UUA', 'CUU', 'UUA', 'CUU'),
 ('UUG', 'CUA', 'UUA', 'CUU', 'UUA', 'CUG'),
 ('UUG', 'CUA', 'UUA', 'CUU', 'UUA', 'CUC'),
 ('UUG', 'CUA', 'UUA', 'CUU', 'UUG', 'CUA'),
 ('UUG', 'CUA', 'UUA', 'CUU', 'UUG', 'CUU'),
 ('UUG', 'CUA', 'UUA', 'CUU', 'UUG', 'CUG'),
 ('UUG', 'CUA', 'UUA', 'CUU', 'UUG', 'CUC'),
 ('UUG', 'CUA', 'UUA', 'CUG', 'UUA', 'CUA'),
 ('UUG', 'CUA', 'UUA', 'CUG', 'UUA', 'CUU'),
 ('UUG', 'CUA', 'UUA', 'CUG', 'UUA', 'CUG'),
 ('UUG', 'CUA', 'UUA', 'CUG', 'UUA', 'CUC'),
 ('UUG', 'CUA', 'UUA', 'CUG', 'UUG', 'CUA'),
 ('UUG', 'CUA', 'UUA', 'CUG', 'UUG', 'CUU'),
 ('UUG', 'CUA', 'UUA', 'CUG', 'UUG', 'CUG'),
 ('UUG', 'CUA', 'UUA', 'CUG', 'UUG', 'CUC'),
 ('UUG', 'CUA', 'UUA', 'CUC', 'UUA', 'CUA'),
 ('UUG', 'CUA', 'UUA', 'CUC', 'UUA', 'CUU'),
 ('UUG', 'CUA', 'UUA', 'CUC', 'UUA', 'CUG'),
 ('UUG', 'CUA', 'UUA', 'CUC', 'UUA', 'CUC'),
 ('UUG', 'CUA', 'UUA', 'CUC', 'UUG', 'CUA'),
 ('UUG', 'CUA', 'UUA', 'CUC', 'UUG', 'CUU'),
 ('UUG', 'CUA', 'UUA', 'CUC', 'UUG', 'CUG'),
 ('UUG', 'CUA', 'UUA', 'CUC', 'UUG', 'CUC'),
 ('UUG', 'CUA', 'UUG', 'CUA', 'UUA', 'CUA'),
 ('UUG', 'CUA', 'UUG', 'CUA', 'UUA', 'CUU'),
 ('UUG', 'CUA', 'UUG', 'CUA', 'UUA', 'CUG'),
 ('UUG', 'CUA', 'UUG', 'CUA', 'UUA', 'CUC'),
 ('UUG', 'CUA', 'UUG', 'CUA', 'UUG', 'CUA'),
 ('UUG', 'CUA', 'UUG', 'CUA', 'UUG', 'CUU'),
 ('UUG', 'CUA', 'UUG', 'CUA', 'UUG', 'CUG'),
 ('UUG', 'CUA', 'UUG', 'CUA', 'UUG', 'CUC'),
 ('UUG', 'CUA', 'UUG', 'CUU', 'UUA', 'CUA'),
 ('UUG', 'CUA', 'UUG', 'CUU', 'UUA', 'CUU'),
 ('UUG', 'CUA', 'UUG', 'CUU', 'UUA', 'CUG'),
 ('UUG', 'CUA', 'UUG', 'CUU', 'UUA', 'CUC'),
 ('UUG', 'CUA', 'UUG', 'CUU', 'UUG', 'CUA'),
 ('UUG', 'CUA', 'UUG', 'CUU', 'UUG', 'CUU'),
 ('UUG', 'CUA', 'UUG', 'CUU', 'UUG', 'CUG'),
 ('UUG', 'CUA', 'UUG', 'CUU', 'UUG', 'CUC'),
 ('UUG', 'CUA', 'UUG', 'CUG', 'UUA', 'CUA'),
 ('UUG', 'CUA', 'UUG', 'CUG', 'UUA', 'CUU'),
 ('UUG', 'CUA', 'UUG', 'CUG', 'UUA', 'CUG'),
 ('UUG', 'CUA', 'UUG', 'CUG', 'UUA', 'CUC'),
 ('UUG', 'CUA', 'UUG', 'CUG', 'UUG', 'CUA'),
 ('UUG', 'CUA', 'UUG', 'CUG', 'UUG', 'CUU'),
 ('UUG', 'CUA', 'UUG', 'CUG', 'UUG', 'CUG'),
 ('UUG', 'CUA', 'UUG', 'CUG', 'UUG', 'CUC'),
 ('UUG', 'CUA', 'UUG', 'CUC', 'UUA', 'CUA'),
 ('UUG', 'CUA', 'UUG', 'CUC', 'UUA', 'CUU'),
 ('UUG', 'CUA', 'UUG', 'CUC', 'UUA', 'CUG'),
 ('UUG', 'CUA', 'UUG', 'CUC', 'UUA', 'CUC'),
 ('UUG', 'CUA', 'UUG', 'CUC', 'UUG', 'CUA'),
 ('UUG', 'CUA', 'UUG', 'CUC', 'UUG', 'CUU'),
 ('UUG', 'CUA', 'UUG', 'CUC', 'UUG', 'CUG'),
 ('UUG', 'CUA', 'UUG', 'CUC', 'UUG', 'CUC'),
 ('UUG', 'CUU', 'UUA', 'CUA', 'UUA', 'CUA'),
 ('UUG', 'CUU', 'UUA', 'CUA', 'UUA', 'CUU'),
 ('UUG', 'CUU', 'UUA', 'CUA', 'UUA', 'CUG'),
 ('UUG', 'CUU', 'UUA', 'CUA', 'UUA', 'CUC'),
 ('UUG', 'CUU', 'UUA', 'CUA', 'UUG', 'CUA'),
 ('UUG', 'CUU', 'UUA', 'CUA', 'UUG', 'CUU'),
 ('UUG', 'CUU', 'UUA', 'CUA', 'UUG', 'CUG'),
 ('UUG', 'CUU', 'UUA', 'CUA', 'UUG', 'CUC'),
 ('UUG', 'CUU', 'UUA', 'CUU', 'UUA', 'CUA'),
 ('UUG', 'CUU', 'UUA', 'CUU', 'UUA', 'CUU'),
 ('UUG', 'CUU', 'UUA', 'CUU', 'UUA', 'CUG'),
 ('UUG', 'CUU', 'UUA', 'CUU', 'UUA', 'CUC'),
 ('UUG', 'CUU', 'UUA', 'CUU', 'UUG', 'CUA'),
 ('UUG', 'CUU', 'UUA', 'CUU', 'UUG', 'CUU'),
 ('UUG', 'CUU', 'UUA', 'CUU', 'UUG', 'CUG'),
 ('UUG', 'CUU', 'UUA', 'CUU', 'UUG', 'CUC'),
 ('UUG', 'CUU', 'UUA', 'CUG', 'UUA', 'CUA'),
 ('UUG', 'CUU', 'UUA', 'CUG', 'UUA', 'CUU'),
 ('UUG', 'CUU', 'UUA', 'CUG', 'UUA', 'CUG'),
 ('UUG', 'CUU', 'UUA', 'CUG', 'UUA', 'CUC'),
 ('UUG', 'CUU', 'UUA', 'CUG', 'UUG', 'CUA'),
 ('UUG', 'CUU', 'UUA', 'CUG', 'UUG', 'CUU'),
 ('UUG', 'CUU', 'UUA', 'CUG', 'UUG', 'CUG'),
 ('UUG', 'CUU', 'UUA', 'CUG', 'UUG', 'CUC'),
 ('UUG', 'CUU', 'UUA', 'CUC', 'UUA', 'CUA'),
 ('UUG', 'CUU', 'UUA', 'CUC', 'UUA', 'CUU'),
 ('UUG', 'CUU', 'UUA', 'CUC', 'UUA', 'CUG'),
 ('UUG', 'CUU', 'UUA', 'CUC', 'UUA', 'CUC'),
 ('UUG', 'CUU', 'UUA', 'CUC', 'UUG', 'CUA'),
 ('UUG', 'CUU', 'UUA', 'CUC', 'UUG', 'CUU'),
 ('UUG', 'CUU', 'UUA', 'CUC', 'UUG', 'CUG'),
 ('UUG', 'CUU', 'UUA', 'CUC', 'UUG', 'CUC'),
 ('UUG', 'CUU', 'UUG', 'CUA', 'UUA', 'CUA'),
 ('UUG', 'CUU', 'UUG', 'CUA', 'UUA', 'CUU'),
 ('UUG', 'CUU', 'UUG', 'CUA', 'UUA', 'CUG'),
 ('UUG', 'CUU', 'UUG', 'CUA', 'UUA', 'CUC'),
 ('UUG', 'CUU', 'UUG', 'CUA', 'UUG', 'CUA'),
 ('UUG', 'CUU', 'UUG', 'CUA', 'UUG', 'CUU'),
 ('UUG', 'CUU', 'UUG', 'CUA', 'UUG', 'CUG'),
 ('UUG', 'CUU', 'UUG', 'CUA', 'UUG', 'CUC'),
 ('UUG', 'CUU', 'UUG', 'CUU', 'UUA', 'CUA'),
 ('UUG', 'CUU', 'UUG', 'CUU', 'UUA', 'CUU'),
 ('UUG', 'CUU', 'UUG', 'CUU', 'UUA', 'CUG'),
 ('UUG', 'CUU', 'UUG', 'CUU', 'UUA', 'CUC'),
 ('UUG', 'CUU', 'UUG', 'CUU', 'UUG', 'CUA'),
 ('UUG', 'CUU', 'UUG', 'CUU', 'UUG', 'CUU'),
 ('UUG', 'CUU', 'UUG', 'CUU', 'UUG', 'CUG'),
 ('UUG', 'CUU', 'UUG', 'CUU', 'UUG', 'CUC'),
 ('UUG', 'CUU', 'UUG', 'CUG', 'UUA', 'CUA'),
 ('UUG', 'CUU', 'UUG', 'CUG', 'UUA', 'CUU'),
 ('UUG', 'CUU', 'UUG', 'CUG', 'UUA', 'CUG'),
 ('UUG', 'CUU', 'UUG', 'CUG', 'UUA', 'CUC'),
 ('UUG', 'CUU', 'UUG', 'CUG', 'UUG', 'CUA'),
 ('UUG', 'CUU', 'UUG', 'CUG', 'UUG', 'CUU'),
 ('UUG', 'CUU', 'UUG', 'CUG', 'UUG', 'CUG'),
 ('UUG', 'CUU', 'UUG', 'CUG', 'UUG', 'CUC'),
 ('UUG', 'CUU', 'UUG', 'CUC', 'UUA', 'CUA'),
 ('UUG', 'CUU', 'UUG', 'CUC', 'UUA', 'CUU'),
 ('UUG', 'CUU', 'UUG', 'CUC', 'UUA', 'CUG'),
 ('UUG', 'CUU', 'UUG', 'CUC', 'UUA', 'CUC'),
 ('UUG', 'CUU', 'UUG', 'CUC', 'UUG', 'CUA'),
 ('UUG', 'CUU', 'UUG', 'CUC', 'UUG', 'CUU'),
 ('UUG', 'CUU', 'UUG', 'CUC', 'UUG', 'CUG'),
 ('UUG', 'CUU', 'UUG', 'CUC', 'UUG', 'CUC'),
 ('UUG', 'CUG', 'UUA', 'CUA', 'UUA', 'CUA'),
 ('UUG', 'CUG', 'UUA', 'CUA', 'UUA', 'CUU'),
 ('UUG', 'CUG', 'UUA', 'CUA', 'UUA', 'CUG'),
 ('UUG', 'CUG', 'UUA', 'CUA', 'UUA', 'CUC'),
 ('UUG', 'CUG', 'UUA', 'CUA', 'UUG', 'CUA'),
 ('UUG', 'CUG', 'UUA', 'CUA', 'UUG', 'CUU'),
 ('UUG', 'CUG', 'UUA', 'CUA', 'UUG', 'CUG'),
 ('UUG', 'CUG', 'UUA', 'CUA', 'UUG', 'CUC'),
 ('UUG', 'CUG', 'UUA', 'CUU', 'UUA', 'CUA'),
 ('UUG', 'CUG', 'UUA', 'CUU', 'UUA', 'CUU'),
 ('UUG', 'CUG', 'UUA', 'CUU', 'UUA', 'CUG'),
 ('UUG', 'CUG', 'UUA', 'CUU', 'UUA', 'CUC'),
 ('UUG', 'CUG', 'UUA', 'CUU', 'UUG', 'CUA'),
 ('UUG', 'CUG', 'UUA', 'CUU', 'UUG', 'CUU'),
 ('UUG', 'CUG', 'UUA', 'CUU', 'UUG', 'CUG'),
 ('UUG', 'CUG', 'UUA', 'CUU', 'UUG', 'CUC'),
 ('UUG', 'CUG', 'UUA', 'CUG', 'UUA', 'CUA'),
 ('UUG', 'CUG', 'UUA', 'CUG', 'UUA', 'CUU'),
 ('UUG', 'CUG', 'UUA', 'CUG', 'UUA', 'CUG'),
 ('UUG', 'CUG', 'UUA', 'CUG', 'UUA', 'CUC'),
 ('UUG', 'CUG', 'UUA', 'CUG', 'UUG', 'CUA'),
 ('UUG', 'CUG', 'UUA', 'CUG', 'UUG', 'CUU'),
 ('UUG', 'CUG', 'UUA', 'CUG', 'UUG', 'CUG'),
 ('UUG', 'CUG', 'UUA', 'CUG', 'UUG', 'CUC'),
 ('UUG', 'CUG', 'UUA', 'CUC', 'UUA', 'CUA'),
 ('UUG', 'CUG', 'UUA', 'CUC', 'UUA', 'CUU'),
 ('UUG', 'CUG', 'UUA', 'CUC', 'UUA', 'CUG'),
 ('UUG', 'CUG', 'UUA', 'CUC', 'UUA', 'CUC'),
 ('UUG', 'CUG', 'UUA', 'CUC', 'UUG', 'CUA'),
 ('UUG', 'CUG', 'UUA', 'CUC', 'UUG', 'CUU'),
 ('UUG', 'CUG', 'UUA', 'CUC', 'UUG', 'CUG'),
 ('UUG', 'CUG', 'UUA', 'CUC', 'UUG', 'CUC'),
 ('UUG', 'CUG', 'UUG', 'CUA', 'UUA', 'CUA'),
 ('UUG', 'CUG', 'UUG', 'CUA', 'UUA', 'CUU'),
 ('UUG', 'CUG', 'UUG', 'CUA', 'UUA', 'CUG'),
 ('UUG', 'CUG', 'UUG', 'CUA', 'UUA', 'CUC'),
 ('UUG', 'CUG', 'UUG', 'CUA', 'UUG', 'CUA'),
 ('UUG', 'CUG', 'UUG', 'CUA', 'UUG', 'CUU'),
 ('UUG', 'CUG', 'UUG', 'CUA', 'UUG', 'CUG'),
 ('UUG', 'CUG', 'UUG', 'CUA', 'UUG', 'CUC'),
 ('UUG', 'CUG', 'UUG', 'CUU', 'UUA', 'CUA'),
 ('UUG', 'CUG', 'UUG', 'CUU', 'UUA', 'CUU'),
 ('UUG', 'CUG', 'UUG', 'CUU', 'UUA', 'CUG'),
 ('UUG', 'CUG', 'UUG', 'CUU', 'UUA', 'CUC'),
 ('UUG', 'CUG', 'UUG', 'CUU', 'UUG', 'CUA'),
 ('UUG', 'CUG', 'UUG', 'CUU', 'UUG', 'CUU'),
 ('UUG', 'CUG', 'UUG', 'CUU', 'UUG', 'CUG'),
 ('UUG', 'CUG', 'UUG', 'CUU', 'UUG', 'CUC'),
 ('UUG', 'CUG', 'UUG', 'CUG', 'UUA', 'CUA'),
 ('UUG', 'CUG', 'UUG', 'CUG', 'UUA', 'CUU'),
 ('UUG', 'CUG', 'UUG', 'CUG', 'UUA', 'CUG'),
 ('UUG', 'CUG', 'UUG', 'CUG', 'UUA', 'CUC'),
 ('UUG', 'CUG', 'UUG', 'CUG', 'UUG', 'CUA'),
 ('UUG', 'CUG', 'UUG', 'CUG', 'UUG', 'CUU'),
 ('UUG', 'CUG', 'UUG', 'CUG', 'UUG', 'CUG'),
 ('UUG', 'CUG', 'UUG', 'CUG', 'UUG', 'CUC'),
 ('UUG', 'CUG', 'UUG', 'CUC', 'UUA', 'CUA'),
 ('UUG', 'CUG', 'UUG', 'CUC', 'UUA', 'CUU'),
 ('UUG', 'CUG', 'UUG', 'CUC', 'UUA', 'CUG'),
 ('UUG', 'CUG', 'UUG', 'CUC', 'UUA', 'CUC'),
 ('UUG', 'CUG', 'UUG', 'CUC', 'UUG', 'CUA'),
 ('UUG', 'CUG', 'UUG', 'CUC', 'UUG', 'CUU'),
 ('UUG', 'CUG', 'UUG', 'CUC', 'UUG', 'CUG'),
 ('UUG', 'CUG', 'UUG', 'CUC', 'UUG', 'CUC'),
 ('UUG', 'CUC', 'UUA', 'CUA', 'UUA', 'CUA'),
 ('UUG', 'CUC', 'UUA', 'CUA', 'UUA', 'CUU'),
 ('UUG', 'CUC', 'UUA', 'CUA', 'UUA', 'CUG'),
 ('UUG', 'CUC', 'UUA', 'CUA', 'UUA', 'CUC'),
 ('UUG', 'CUC', 'UUA', 'CUA', 'UUG', 'CUA'),
 ('UUG', 'CUC', 'UUA', 'CUA', 'UUG', 'CUU'),
 ('UUG', 'CUC', 'UUA', 'CUA', 'UUG', 'CUG'),
 ('UUG', 'CUC', 'UUA', 'CUA', 'UUG', 'CUC'),
 ('UUG', 'CUC', 'UUA', 'CUU', 'UUA', 'CUA'),
 ('UUG', 'CUC', 'UUA', 'CUU', 'UUA', 'CUU'),
 ('UUG', 'CUC', 'UUA', 'CUU', 'UUA', 'CUG'),
 ('UUG', 'CUC', 'UUA', 'CUU', 'UUA', 'CUC'),
 ('UUG', 'CUC', 'UUA', 'CUU', 'UUG', 'CUA'),
 ('UUG', 'CUC', 'UUA', 'CUU', 'UUG', 'CUU'),
 ('UUG', 'CUC', 'UUA', 'CUU', 'UUG', 'CUG'),
 ('UUG', 'CUC', 'UUA', 'CUU', 'UUG', 'CUC'),
 ('UUG', 'CUC', 'UUA', 'CUG', 'UUA', 'CUA'),
 ('UUG', 'CUC', 'UUA', 'CUG', 'UUA', 'CUU'),
 ('UUG', 'CUC', 'UUA', 'CUG', 'UUA', 'CUG'),
 ('UUG', 'CUC', 'UUA', 'CUG', 'UUA', 'CUC'),
 ('UUG', 'CUC', 'UUA', 'CUG', 'UUG', 'CUA'),
 ('UUG', 'CUC', 'UUA', 'CUG', 'UUG', 'CUU'),
 ('UUG', 'CUC', 'UUA', 'CUG', 'UUG', 'CUG'),
 ('UUG', 'CUC', 'UUA', 'CUG', 'UUG', 'CUC'),
 ('UUG', 'CUC', 'UUA', 'CUC', 'UUA', 'CUA'),
 ('UUG', 'CUC', 'UUA', 'CUC', 'UUA', 'CUU'),
 ('UUG', 'CUC', 'UUA', 'CUC', 'UUA', 'CUG'),
 ('UUG', 'CUC', 'UUA', 'CUC', 'UUA', 'CUC'),
 ('UUG', 'CUC', 'UUA', 'CUC', 'UUG', 'CUA'),
 ('UUG', 'CUC', 'UUA', 'CUC', 'UUG', 'CUU'),
 ('UUG', 'CUC', 'UUA', 'CUC', 'UUG', 'CUG'),
 ('UUG', 'CUC', 'UUA', 'CUC', 'UUG', 'CUC'),
 ('UUG', 'CUC', 'UUG', 'CUA', 'UUA', 'CUA'),
 ('UUG', 'CUC', 'UUG', 'CUA', 'UUA', 'CUU'),
 ('UUG', 'CUC', 'UUG', 'CUA', 'UUA', 'CUG'),
 ('UUG', 'CUC', 'UUG', 'CUA', 'UUA', 'CUC'),
 ('UUG', 'CUC', 'UUG', 'CUA', 'UUG', 'CUA'),
 ('UUG', 'CUC', 'UUG', 'CUA', 'UUG', 'CUU'),
 ('UUG', 'CUC', 'UUG', 'CUA', 'UUG', 'CUG'),
 ('UUG', 'CUC', 'UUG', 'CUA', 'UUG', 'CUC'),
 ('UUG', 'CUC', 'UUG', 'CUU', 'UUA', 'CUA'),
 ('UUG', 'CUC', 'UUG', 'CUU', 'UUA', 'CUU'),
 ('UUG', 'CUC', 'UUG', 'CUU', 'UUA', 'CUG'),
 ('UUG', 'CUC', 'UUG', 'CUU', 'UUA', 'CUC'),
 ('UUG', 'CUC', 'UUG', 'CUU', 'UUG', 'CUA'),
 ('UUG', 'CUC', 'UUG', 'CUU', 'UUG', 'CUU'),
 ('UUG', 'CUC', 'UUG', 'CUU', 'UUG', 'CUG'),
 ('UUG', 'CUC', 'UUG', 'CUU', 'UUG', 'CUC'),
 ('UUG', 'CUC', 'UUG', 'CUG', 'UUA', 'CUA'),
 ('UUG', 'CUC', 'UUG', 'CUG', 'UUA', 'CUU'),
 ('UUG', 'CUC', 'UUG', 'CUG', 'UUA', 'CUG'),
 ('UUG', 'CUC', 'UUG', 'CUG', 'UUA', 'CUC'),
 ('UUG', 'CUC', 'UUG', 'CUG', 'UUG', 'CUA'),
 ('UUG', 'CUC', 'UUG', 'CUG', 'UUG', 'CUU'),
 ('UUG', 'CUC', 'UUG', 'CUG', 'UUG', 'CUG'),
 ('UUG', 'CUC', 'UUG', 'CUG', 'UUG', 'CUC'),
 ('UUG', 'CUC', 'UUG', 'CUC', 'UUA', 'CUA'),
 ('UUG', 'CUC', 'UUG', 'CUC', 'UUA', 'CUU'),
 ('UUG', 'CUC', 'UUG', 'CUC', 'UUA', 'CUG'),
 ('UUG', 'CUC', 'UUG', 'CUC', 'UUA', 'CUC'),
 ('UUG', 'CUC', 'UUG', 'CUC', 'UUG', 'CUA'),
 ('UUG', 'CUC', 'UUG', 'CUC', 'UUG', 'CUU'),
 ('UUG', 'CUC', 'UUG', 'CUC', 'UUG', 'CUG'),
 ('UUG', 'CUC', 'UUG', 'CUC', 'UUG', 'CUC')]
len(x)
512/3
len("LPEETMDFMVQHTSFKEMKKNPMTNYTTVPQELMDHSISPFMRKGM")
46 * 3
LGHYUJUSOSKKO
list1 = []
def dna_num (rna_seq):
    for i in range(0,len(rna_seq)-2):
        list1.append(rna_seq[i:i+3])
    return list1, len(list1)

a = "LPEETMDFMVQHTSFKEMKKNPMTNYTTVPQELMDHSISPFMRKGM"
a[0:-2]

dna_num ("LPEETMDFMVQHTSFKEMKKNPMTNYTTVPQELMDHSISPFMRKGM")
