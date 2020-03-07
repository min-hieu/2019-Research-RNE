import eel
from Super_Func import *
from Super_Seq import *

eel.init('web')

@eel.expose
def testing(seq,site,thres):

    thres = int(thres)

    #20-23+27

    # MLPGVGLTPS AAQTARQHPK MHLAHSTLKP AAHLIGDPSK QNSLLWRANT
    # DRAFLQDGFS LSNNSLLVPT SGIYFVYSQV VFSGKAYSPK ATSSPLYLAH
    # EVQLFSSQYP FHVPLLSSQK MVYPGLQEPW LHSMYHGAAF QLTQGDQLST
    # HTDGIPHLVL SPSTVFFGAF AL

    protein = Super_Seq(seq, name="TNFB")
    site = binding_site_converter(site)

    charged_protein = supercharge_arbitrary(protein, threshold=thres, charge=1, binding_site=site)

    protein.compare_lcd(charged_protein, threshold=thres)

    return charged_protein.seq
    

my_options = {
    'mode': "edge", #or "chrome-app",
    'size': (550, 650),
    'port': 8080,
    'chromeFlags': ["--start-fullscreen", "--browser-startup-dialog"]
}

eel.start('index.html', size=(550, 650), port=8080)
while True:
    eel.sleep(10)