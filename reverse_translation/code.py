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
