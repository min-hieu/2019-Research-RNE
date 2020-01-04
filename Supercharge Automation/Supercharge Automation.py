aa = [
    'D','T','S','E','P','G','A','C','V','M','I','L','Y','F','H','K','R','W','Q','N'
]
f = {
    'D': -1, 'T': 0, 'S': 0, 'E': -1, 'P': 0, 'G': 0, 'A': 0, 'C': -0.1, 'V': 0, 'M': 0, 'I': 0, 'L': 0, 'Y': 0,
    'F': 0, 'H': 0.1, 'K': 1, 'R': 1, 'W': 0, 'Q': 0, 'N': 0
}

def delSpace(seq):
    n = len(seq)
    newStr = ""
    aa = ['D','T','S','E','P','G','A','C','V','M','I','L','Y','F','H','K','R','W','Q','N']
    for i in range(n):
        if seq[i] in aa:
            newStr += seq[i]
    return newStr


def scoreCharge(seq, size=20):

    seq = delSpace(seq)
    w = size

    indivisualCharge = []

    for i in range(len(seq) - (w - 1)):
        sum_ = 0
        for j in range(i, i + w):
            sum_ += f[seq[j]]
        indivisualCharge.append(sum_)

    return indivisualCharge

def scoreConsurf(filename):
    with open(filename,'r') as f:




if __name__ == "__main__":
    intp = """
    MLPGVGLTPSAAQTARQHPKMHLAHSTLKPAAHLIGDPSKQNSLLWRANTDRAFLQDGFSLSNNSLLVPTSGIYFVY
    SQVVFSGKAYSPKATSSPLYLAHEVQLFSSQYPFHVPLLSSQKMVYPGLQEPWLHSMYHGAAFQLTQGDQLSTHTDG
    IPHLVLSPSTVFFGAFAL      
    """
    scoreCharge(intp)