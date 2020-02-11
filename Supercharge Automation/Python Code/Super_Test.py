
from Super_Func import *
from Super_Seq import *

def console():
    # seq1 = input("1st protein:  ")
    # name1 = input("1st protein name:  ")
    # binding_site = input("bindind site: ")

    seq_sample = """MLPGVGLTPS AAQTARQHPK MHLAHSTLKP AAHLIGDPSK QNSLLWRANT
    DRAFLQDGFS LSNNSLLVPT SGIYFVYSQV VFSGKAYSPK ATSSPLYLAH
    EVQLFSSQYP FHVPLLSSQK MVYPGLQEPW LHSMYHGAAF QLTQGDQLST
    HTDGIPHLVL SPSTVFFGAF AL"""

    threshold = 1
    
    protein = Super_Seq(seq_test, name="TNFB")
    binding_site = binding_site_converter("20-23+27")
    charged_protein = supercharge_arbitrary(protein, threshold=threshold, charge=1, binding_site=binding_site)

    protein.compare_lcd(charged_protein, threshold=threshold)


console()


