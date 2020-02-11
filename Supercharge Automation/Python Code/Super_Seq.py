
import math
import matplotlib.pyplot as plt


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
            charges.append(round(sum_,3))

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
        plt.plot(x, y2, color='orange')
        plt.axhline(threshold, 0, 1, alpha=0.6, color="#00F481")
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

