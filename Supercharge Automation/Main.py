"""
Author: Nguyen Minh Hieu




"""

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

        score = list()

        with open(filename, "r") as grade:
            for i, line in enumerate(grade):
                if i > 14: score.append(line.split()[4])

        return score

    def __init__(self, sequence, filename):
        self.__raw_seq = sequence
        self.__seq = self.format_seq(sequence)
        self.__charge = self.seq_charge(sequence)
        self.__consurf = self.consurf_reader(filename)

    def consurf(self):

        return None

    def min_lcd(self, threshold = 4):
        seq = self.__seq
        charge = self.__charge
        for i in range(len(charge)):
            if i > threshold:
