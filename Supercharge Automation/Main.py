"""
Author: Nguyen Minh Hieu




"""

import math
import matplotlib.pyplot as plt
import random

class Super_Seq:

    chargeDict = {'D': -1, 'T': 0, 'S': 0, 'E': -1, 'P': 0, 'G': 0, 'A': 0, 'C': -0.1, 'V': 0, 'M': 0, 'I': 0, 'L': 0, 'Y': 0,
                  'F': 0, 'H': 0.1, 'K': 1, 'R': 1, 'W': 0, 'Q': 0, 'N': 0}

    aa = ['D', 'T', 'S', 'E', 'P', 'G', 'A', 'C', 'V', 'M', 'I', 'L', 'Y', 'F', 'H', 'K', 'R', 'W', 'Q', 'N']

    @staticmethod
    def format_seq(sequence):
        newStr = ""
        sequence.upper()
        for i in sequence:
            if i in Super_Seq.aa: newStr += i

        return newStr

    @staticmethod
    def seq_charge(sequence, window_size=20):

        seq = Super_Seq.format_seq(sequence)
        w = window_size
        charges = []

        for i in range(len(seq) - (w - 1)):
            sum_ = 0
            for j in range(i, i + w):
                sum_ += Super_Seq.chargeDict[seq[j]]
            charges.append(sum_)

        return charges

    @staticmethod
    def consurf_reader(filename):
        assert type(filename) == str
        if filename == "": return None

        score = list()

        with open(filename, "r") as grade:
            for i, line in enumerate(grade):
                if i > 14: score.append(line.split()[4])

        return score

    def __init__(self, sequence, filename="", name="protein of interest"):
        self.__raw_seq = sequence
        self.seq = self.format_seq(sequence)
        self.charge = self.seq_charge(sequence)
        self.consurf = self.consurf_reader(filename)
        self.name = name
        self.location = filename

    def consurf(self):
        return None

    def display_lcd(self):
        charge = self.charge

        x = range(len(charge))
        y = charge

        plt.plot(x, y)
        plt.axhline(4,0,1,alpha = 0.6, color = "orange")
        plt.xlabel("Window Index")
        plt.ylabel("Charge Index")
        plt.suptitle(self.name, fontsize = 20)

        plt.show()

    def compare_lcd(self, protein, threshold=4):
        assert type(protein) == Super_Seq
        x = range(len(self.charge))
        y1 = self.charge
        y2 = protein.charge

        title = self.name + " normal(blue) vs charged(orange)"

        plt.plot(x, y1, color='blue')
        plt.plot(x, y2, color='red')
        plt.axhline(threshold, 0, 1, alpha=0.6, color="black")
        plt.xlabel("Window Index")
        plt.ylabel("Charge Index")
        plt.suptitle(title, fontsize=20)

        plt.show()

    def comp_sequence(self, seq2):
        square = 0
        n = len(self.seq)

        for i in range(n):
            square += (self.chargeDict[self.seq[i]] - self.chargeDict[seq2.seq[i]])**2

        return math.sqrt(sum / n)

    # def min_lcd(self, threshold = 4):
    #     seq = self.seq
    #     charge = self.charge
    #     for i in range(len(charge)):
    #         if i > threshold:


def binding_site_converter(site):
    assert type(site) == str
    assert len(site) > 0
    # format: #1-#2+#3 (from #1 to #2 add #3)
    # output (tuple): (#1,...,#2,#3)

    siteList = site.split("+")
    extendedList = []
    to_remove = []
    for elem in siteList:
        if "-" in elem:
            start, end = map(int, elem.split("-"))
            to_remove.append(elem)
            extendedList += list(range(start, end+1))

    for elem in to_remove:
        siteList.remove(elem)

    siteList = list(map(int, siteList)) + extendedList

    return tuple(siteList)


def superchargeMethod(seq, method_tuple,threshold, charge):
    #output a tuple of indexes based on the given method
    assert type(seq) == Super_Seq
    assert sum(method_tuple) > 0

    LCD, consf = method_tuple
    # if LCD and consf:
    #
    # elif LCD:
    # elif consf:
    # else:


def supercharge_arbitrary(seq, threshold, charge, binding_site):
    assert type(seq) == Super_Seq
    sequence = list(seq.seq)
    ch = seq.charge

    ok = len(ch)
    ok2 = len(sequence)

    for i in range(len(ch)):
        if ch[i]*charge > threshold: 
            for j in range(20):
                if sequence[i+j] in ["R", "K"] and i+j not in binding_site:
                    sequence[i+j] = random.choice(["D", "E"])
                    curr_sequence = "".join(sequence[i:i+21])
                    curr_charge = Super_Seq.seq_charge(curr_sequence)[0]
                    if curr_charge <= threshold: break

    sequence = ''.join(sequence)
                    
    return Super_Seq(sequence, seq.location, seq.name)


def supercharge(protein, binding_site, method, threshold, charge):
    # assert [type(protein), type(binding_site), type(method), type(threshold), type(charge)] \
    #        == [Super_Seq, str, tuple, float, str]

    binding_site = binding_site_converter(binding_site)
    charged_protein = supercharge_arbitrary(protein, threshold, charge, binding_site)
    protein.compare_lcd(charged_protein)


# seq_test = """MLPGVGLTPS AAQTARQHPK MHLAHSTLKP AAHLIGDPSK QNSLLWRANT
# DRAFLQDGFS LSNNSLLVPT SGIYFVYSQV VFSGKAYSPK ATSSPLYLAH
# EVQLFSSQYP FHVPLLSSQK MVYPGLQEPW LHSMYHGAAF QLTQGDQLST
# HTDGIPHLVL SPSTVFFGAF AL"""
#
# TNFB = Super_Seq(seq_test, "", "TNFB")
# TNFB.display_lcd()


def console():
    # seq1 = input("1st protein:  ")
    # name1 = input("1st protein name:  ")
    # binding_site = input("bindind site: ")

    seq_test = """MLPGVGLTPS AAQTARQHPK MHLAHSTLKP AAHLIGDPSK QNSLLWRANT
    DRAFLQDGFS LSNNSLLVPT SGIYFVYSQV VFSGKAYSPK ATSSPLYLAH
    EVQLFSSQYP FHVPLLSSQK MVYPGLQEPW LHSMYHGAAF QLTQGDQLST
    HTDGIPHLVL SPSTVFFGAF AL"""

    threshold = 1
    
    protein = Super_Seq(seq_test, name="TNFB")
    binding_site = binding_site_converter("20-23+27")
    charged_protein = supercharge_arbitrary(protein, threshold=threshold, charge=1, binding_site=binding_site)

    protein.compare_lcd(charged_protein, threshold=threshold)


console()
