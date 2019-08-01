# Python 3.7
# encoded in UTF-8
# author Hyunjong Byun
# Version 1.0:
# Does the same thing as Choi-Cho's PrimerBuilder, but written in python and can handle multiple sequences in FASTA format.
# Version 1.1:
# Tagging functionality added.
# Version 1.2:
# GUI support using tkinter added.

from SequenceProcessor import DNAset, DNAXset, RNAset, RNAXset
from SequenceProcessor import fasta_reader, fasta_writer
from SequenceProcessor import inverse_dictionary, mode_set
from SequenceProcessor import reverse_complement

enzset = { 'XbaI':'TCTAGA', 'NdeI':'CATATG', 'KpnI':'GGTACC', 'SacI':'GAGCTC', \
           'HindIII':'AAGCTT', 'EcoRI':'GAATTC', 'SacII':'CCGCGG' }

tagset = { 'His6':'CATCACCATCACCATCAC', 'His6_alt':'CACCACCACCACCACCAC' }



def primer_builder(fr, item, enz, tag, Tm=60, mode='DNA'):
    """
    takes a specifier for forward or reverse: 'f' for forward, 'r' for reverse
    takes a single tuple of (header, sequence) as the input.
    also takes the name of forward or reverse primer.
    returns list of tuples (header, sequence)
    naming convention of AhnLab is used to name the primers.
    """
    if not enz == 'None':
        if not item[1].find(enzset[enz]) == -1:
            print('Error: sequence ' + item[0] + ' contains the restriction site of ' + enz + '.')
            return
    temp_header = item[0]
    if enz == 'None':
        if fr == 'f':
            temp_header = 'f_' + temp_header
        else:
            temp_header = 'r_' + temp_header
    else:
        if fr == 'f':
            temp_header = enz + '-' + temp_header
        else:
            temp_header = temp_header + '-' + enz
    
    if tag == 'None':
        temp_tag = ''
    elif tag in tagset:
        temp_tag = tagset[tag]
    else:
        temp_tag = tag
    if fr == 'r' and tag != 'None':
        temp_tag = reverse_complement(temp_tag)
    
    calc_Tm = 0
    temp_seq = ''
    if fr == 'f':
        for i in range(len(item[1])):
            if item[1][i] in ('A', 'T', 'U', 'W'):
                calc_Tm += 2
            elif item[1][i] in ('G', 'C', 'S'):
                calc_Tm += 4
            elif item[1][i] in ('R', 'Y', 'N'):
                calc_Tm += 3
            elif item[1][i] in ('B', 'V'):
                calc_Tm += 3.3333
            elif item[1][i] in ('D', 'H'):
                calc_Tm += 2.6667
            temp_seq = temp_seq + item[1][i]
            if calc_Tm >= Tm:
                break
        if not enz == 'None':
            temp_seq = 'GGG' + enzset[enz] + temp_tag + temp_seq
    else:
        for i in range(len(item[1])):
            index = len(item[1])-1-i
            if item[1][index] in ('A', 'T', 'U', 'W'):
                calc_Tm += 2
            elif item[1][index] in ('G', 'C', 'S'):
                calc_Tm += 4
            elif item[1][index] in ('R', 'Y', 'N'):
                calc_Tm += 3
            elif item[1][index] in ('B', 'V'):
                calc_Tm += 3.3333
            elif item[1][index] in ('D', 'H'):
                calc_Tm += 2.6667
            temp_seq = item[1][index] + temp_seq
            if calc_Tm >= Tm:
                break
        temp_seq = reverse_complement(temp_seq)
        if not enz == 'None':
            temp_seq = 'GGG' + reverse_complement(enzset[enz]) + temp_tag + temp_seq
    return temp_header, temp_seq


def _mode_select():
    mode = input('Input mode: ').strip()
    if not mode in ('DNA','DNAX'):
        print('Error: Invalid mode.')
        mode = _mode_select()
    return mode


def _Tm_select():
    Tm = input('Input Tm: ').strip()
    if not Tm.isnumeric():
        return 'Error: Invalid Tm.'
        _Tm_select()
    Tm = int(Tm)
    '''
    if Tm < 50:
        yn = input('Warning: Too low Tm. Wish to continue? (y/n): ')
        if yn == 'y':
            pass
        else:
            Tm = _Tm_select()
    '''
    return Tm



def _primer_select(fr):
    if fr == 'f':
        text = 'forward'
        Text = 'Forward'
    elif fr == 'r':
        text = 'reverse'
        Text = 'Reverse'
    primer = input('Input ' + text + ' restriction site: ').strip()
    if primer == '' or primer == 'None':
        return 'None'
    elif not primer in enzset:
        print('Error: ' + Text + ' primer does not exist in dictionary "enzset".')
        primer = _primer_select(fr)
    return primer



def _tag_select(mode):
    tag = input('Input tag: ').strip()
    if tag == '' or tag == 'None':
        tag == 'None'
    else:
        if tag in tagset:
            pass #Just a pass. Do not implement anything.
        elif tag.upper() in ('HIS6', 'HIS-TAG', 'HISTAG', 'POLYHISTIDINE'):
            tag = 'His6'
        else:
            print('Warning:invalid tag name...')
            yn = input('Wish to handle the input as sequence literal? (y/n): ')
            if yn == 'y':
                comset = mode_set(mode)
                temp_tag = ''
                for char in tag:
                    if char in comset:
                        temp_tag += char
                    else:
                        print('Warning:Character \'' + char + '\' is not a valid character in current mode and was ignored.')
                tag = temp_tag
            else:
                tag = _tag_select(mode)
    return tag



def _tag_fr_select():
    tagfr = input('Tag on forward or reverse? (f/r): ').strip()
    if tagfr == 'f' or tagfr == 'r':
        return tagfr
    else:
        print('Error: Invalid input.')
        tag_fr = _tag_fr_select()  


'''
def _main():
    window = tk.Tk()
    
    window.title('PrimerBuilder Reborn v1.2')
    window.geometry("640x400+100+100")
    window.resizable(True, True)
    
    label_info = tk.Label(window, text='PrimerBuilder Reborn v1.2 by Hyunjong Byun\npowered by SequenceProcessor v2.2 by Hyunjong Byun', \
                          anchor='ne', justify = 'right')
    label_info.pack(padx=10, )
    button_build = tk.Button(window, overrelief='sunken', width=15, command=_build, repeatdelay=1000, repeatinterval=100)
    button_build.pack()
    
    window.mainloop()
'''


def _build():
    pass



def _console_main():
    mode = _mode_select()
    forward = _primer_select('f')
    reverse = _primer_select('r')
    Tm = _Tm_select()
    tag = _tag_select(mode)
    
    if tag != 'None':
        tagfr = _tag_fr_select()
    else:
        tagfr = 'f'
    
    if tagfr == 'f':
        ftag = tag
        rtag = 'None'
    else:
        rtag = tag
        ftag = 'None'        
    print('Processing data...')
    data = fasta_reader(mode)
    processed = []
    for item in data:
        processed.append(primer_builder('f', item, forward, ftag, Tm, mode))
        processed.append(primer_builder('r', item, reverse, rtag, Tm, mode))
    fasta_writer(processed, 99999999, False)
    print('Task finished.')
    input('Press any key to exit.')
    return



if __name__ == '__main__':
    # execute only if run as a script
    _console_main()
