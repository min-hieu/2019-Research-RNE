"""
Python 3.6.8
encoded in UTF-8
author Nguyen Minh Hieu
Version 0.0.0a
"""

import random
from Super_Seq import *


def get_curr_charge(i, seq_list):
    curr_sequence = "".join(seq_list[i:i + 21])
    curr_charge = Super_Seq.seq_charge(curr_sequence)[0]
    return curr_charge


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


def supercharge_arbitrary(seq, threshold, charge, binding_site):
    assert type(seq) == Super_Seq
    sequence = list(seq.seq)
    ch = seq.charge

    for i in range(len(ch)):
        if ch[i]*charge > threshold: #if negatively charge the sequence, charge < 0
            if i == 0:
                for j in range(20):
                    if sequence[j] in ["R", "K"] and j not in binding_site:
                        sequence[j] = random.choice(["D", "E"])
                        curr_charge = get_curr_charge(0, sequence)
                        if curr_charge <= threshold:
                            ch = Super_Seq.seq_charge("".join(sequence))
                            break
            else:
                if sequence[i+20] in ["R", "K"] and i+20 not in binding_site:
                    sequence[i+20] = random.choice(["D", "E"])


    sequence = ''.join(sequence)
                    
    return Super_Seq(sequence, seq.location, seq.name)


def supercharge(protein, binding_site, method, threshold, charge):
    # assert [type(protein), type(binding_site), type(method), type(threshold), type(charge)] \
    #        == [Super_Seq, str, tuple, float, str]

    binding_site = binding_site_converter(binding_site)
    charged_protein = supercharge_arbitrary(protein, threshold, charge, binding_site)
    protein.compare_lcd(charged_protein)



def supercharge_exhaustive(protein, binding_site, method, threshold, charge):
    assert type(seq) == Super_Seq
    ch = seq.charge
    ch_aa = random.choice(["D","E"])

    for k in range(len(ch)):
        score = 0
        sequence = list(seq.seq)
        for i in range(len(ch)):
            if ch[i]*charge > threshold: #if negatively charge the sequence, charge < 0
                if i == k:
                    for j in range(20):
                        if sequence[j] in ["R", "K"] and j not in binding_site:
                            sequence[j] = ch_aa
                            curr_charge = get_curr_charge(0, sequence)
                            if curr_charge <= threshold:
                                ch = Super_Seq.seq_charge("".join(sequence))
                                break
                else:
                    if sequence[i+20] in ["R", "K"] and i+20 not in binding_site:
                        sequence[i+20] = ch_aa
                        ch = Super_Seq.seq_charge(sequence)

        charged = sequence = ''.join(sequence)
        pass


                    
    return None