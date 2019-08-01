# Python 3.7
# encoded in UTF-8
# author Hyunjong Byun
# Version 2.1:
# Supports for degenerate DNA codes are added. More Object-oriented: one function now does one thing.
# Version 2.2:
# Stop codon now follows the IUPAC notation: '*', not '-'.
# Version 2.3:
# The code pattern is now much more object-oriented, e.g. reverse_complement(). Bottlenecks in fasta_reader() and fasta_writer() are now removed.
import collections

DNAset = {'A', 'T', 'G', 'C'}
DNAXset = {'A', 'T', 'G', 'C', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'}
RNAset = {'A', 'U', 'G', 'C'}
RNAXset = {'A', 'U', 'G', 'C', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'}
AAset = {'G', 'A', 'V', 'L', 'I', 'M', 'P', 'W', 'F', 'S', 'T', 'C', 'Y',\
             'N', 'Q', 'D', 'E', 'K', 'R', 'H', '*'}
AAXset = {'G', 'A', 'V', 'L', 'I', 'M', 'P', 'W', 'F', 'S', 'T', 'C', 'Y',\
              'N', 'Q', 'D', 'E', 'K', 'R', 'H', 'B', 'Z', 'X', '*'}

codonTable = {\
    'TTT':'F', 'TTC':'F',\
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',\
    'ATT':'I', 'ATC':'I', 'ATA':'I',\
    'ATG':'M',\
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',\
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',\
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',\
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',\
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',\
    'TAT':'Y', 'TAC':'Y',\
    'TAA':'*', 'TAG':'*', 'TGA':'*',\
    'CAT':'H', 'CAC':'H',\
    'CAA':'Q', 'CAG':'Q',\
    'AAT':'N', 'AAC':'N',\
    'AAA':'K', 'AAG':'K',\
    'GAT':'D', 'GAC':'D',\
    'GAA':'E', 'GAG':'E',\
    'TGT':'C', 'TGC':'C',\
    'TGG':'W',\
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',\
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

basePairingTable = {\
    'A':'T', 'T':'A', 'G':'C', 'C':'G',\
    'R':'Y', 'Y':'R', 'K':'M', 'M':'K',\
    'S':'S', 'W':'W',\
    'B':'V', 'V':'B', 'D':'H', 'H':'D',\
    'N':'N'}



def inverse_dictionary(f):
    d = {}
    for key in f:
        if f[key] in d:
            d[f[key]].append(key)
        else:
            d[f[key]] = [key]
    return d



def mode_set(mode):
    comset = AAXset #a fail-safe
    if mode == 'DNA':
        comset = DNAset
    elif mode == 'DNAX':
        comset = DNAXset
    elif mode == 'RNA':
        comset = RNAset
    elif mode == 'RNAX':
        comset = RNAXset
    elif mode == 'AA':
        comset = AAset
    elif mode == 'AAX':
        comset = AAXset
    return comset



def fasta_reader(mode, fasta_name='./input.txt'): 
    """
    reads a fasta file. yield list of tuples of strings [(header, sequence),*].
    mode == "DNA" or "DNAX"(DNA with degenerate codes) or "RNA" or "RNAX" or "AA"(protein) or "AAX".
    default of mode is AAX.
    """
    f = open(fasta_name,'r')

    comset = mode_set(mode)
    out = []
    lines = f.readlines()
    started = False
    temp_header = ''
    temp_seq = []
    for line in lines:
        line_lstrip = line.lstrip()
        if line_lstrip.strip() == '':
            continue
        elif line_lstrip[0] == '>':
            if started:
                temp_seq = ''.join(temp_seq_list)
                out.append((temp_header, temp_seq))
            else:
                started = True
            temp_header = line[1:].rstrip()
            temp_seq_list = []
        elif line_lstrip[0] == ';':
            continue
        else:
            l = line_lstrip.strip().upper()
            for i in range(len(l)):
                if l[i] in comset:
                    temp_seq_list.extend(l[i])
                elif l[i].isalpha():
                    continue
                else:
                    print('Error: character \"' + l[i] + '\" does not match with specified type \"'\
                          + mode + '\" in function fasta_reader(' + fasta_name + ').')
    temp_seq = ''.join(temp_seq_list)
    out.append((temp_header, temp_seq))
    f.close()
    return out



def fasta_writer(input_list, length=100, numbered=True, vertically_spaced=True):
    """
    given a list of tuples of strings in form [(header, sequence),*], yields fasta-formatted .txt file.
    line length is given, or can be 100 to be default.
    numbered = True writes numbers in the front only if len(seq) > length.
    numbered = False does not write the numbers in the front.
    vertically_spaced = True writes an additional '\n' between each fasta elements.
    """
    f = open("./output.txt", 'w')
    length = int(length) #a failsafe
    for item in input_list:
        pointer = 0
        f.write('>' + item[0] + '\n')
        numberwriting = False
        if numbered and len(item[1]) > length:
            numberwriting = True
        m = int(len(item[1])/length)
        n = m*length
        for i in range(m + 1):
            if numberwriting:
                if n >= 1000000 and i*length / 1000000 < 1: #trust python's compiler. but ordering might help.
                    f.write(' ')                    
                if n >= 100000 and i*length / 100000 < 1:
                    f.write(' ')                   
                if n >= 10000 and i*length / 10000 < 1:
                    f.write(' ')                
                if n >= 1000 and i*length / 1000 < 1:
                    f.write(' ')
                if n >= 100 and i*length / 100 < 1:
                    f.write(' ')
                if n >= 10 and i*length / 10 < 1:
                    f.write(' ')
                f.write(str(i*length + 1))
                f.write(' ')
            f.write(item[1][i*length:(i+1)*length] + '\n')
        if vertically_spaced:
            f.write('\n')
    f.close()
    return



################## sequence_processor() elements ################
def identity_function(char): #identical to A_to_B_function({}), but this one is faster.
    return char

def A_to_B_function(char_dic): #this function returns a function.
    def result_function(char):
        if char in char_dic:
            return char_dic[char]
        else:
            return char
    return result_function

dic_U_to_T = { 'U':'T' }
dic_T_to_U = { 'T':'U' }
################



def sequence_processor(seq, comset_name, replacer_function=identity_function, not_reverse=True):
    out = collections.deque() #use deque to support extendleft().
    comset = mode_set(comset_name)
    if not_reverse:
        for char in seq:
            if char in comset:
                out.extend(replacer_function(char))
            else:
                print('Error: character \"'+char+'\" in seq in function trascript(seq) is not '+comset_name+' type.')
    else:
        for char in seq:
            if char in comset:
                out.extendleft(replacer_function(char))
            else:
                print('Error: character \"'+char+'\" in seq in function trascript(seq) is not '+comset_name+' type.')        
    return ''.join(out)



def reverse_complement(seq):
    return sequence_processor(seq, 'DNAX', A_to_B_function(basePairingTable), not_reverse = False)

def reverse(seq):
    return sequence_processor(seq, 'AAX', not_reverse = False)

def transcript(seq):
    return sequence_processor(seq,'RNAX', A_to_B_function(dic_U_to_T))

def reverse_transcript(seq):
    return sequence_processor(seq,'DNAX', A_to_B_function(dic_T_to_U))



def translate(data): #Only supports DNA and RNA format. receives list of tuples [(header, sequence),*]
    out = []
    for item in data:
        for readingframe in [0,1,2]:
            temp_translate_list = []
            codon = ''#since max length of codons is only 3, little need to use list.
            for i in range(len(item[1])-readingframe):
                if item[1][i+readingframe] == 'U':
                    codon = codon + 'T'
                else:
                    codon = codon + item[1][i+readingframe]
                if len(codon) == 3:
                    try:
                        temp_translate_list.append(codonTable[codon])
                    except KeyError:
                        print('Error: Non-DNA character in the input in function translate(data).')
                    codon = ''
            temp_translate = ''.join(temp_translate_list)
            out.append((item[0] + '_' + str(readingframe) , temp_translate))
    return out



def translateX(data): #Incomplete. Supports translation of DNAX and RNAX format. Slower than translate().
    out = []
    for item in data:
        for readingframe in [0,1,2]:
            temp_translate = ''
            codon = ''
            for i in range(len(item[1])-readingframe):
                if item[1][i+readingframe] == 'U':
                    codon = codon + 'T'
                elif item[1][i+readingframe] in DNAXset:
                    codon = codon + item[1][i+readingframe]
                if len(codon) == 3:
                    try:
                        temp_translate = temp_translate + codonTable[codon]
                    except KeyError:
                        print('Error: Non-DNA character in the input in function translateX(data).')
                    codon = ''
            out.append((item[0] + '_' + str(readingframe) , temp_translate))
    return out



def mutation_reader(): #Incomplete. 
    pass




def point_mutator(data, mutlist): #Incomplete. 
    """
    input is list of tuples (header, sequence) and list of mutations tuples (position, result)
    """
    pass



def codon_statistics(): #Incomplete. 
    reader("DNA")
    codons = {}
    for i in codonTable.keys():
        codons[i] = 0
    memory = open("./memory.txt","r")
    result = open("./result.csv","w")
    result.write('\"Amino Acid\",\"Codons\",\"Value\"\n')
    line = memory.readline().strip()
    i = 0
    totalAA = 0
    while(i < len(line)-3):
        frame = line[i] + line[i+1] + line[i+2]
        codons[frame] += 1
        totalAA += 1
        i += 3
    codonList = inverse_dictionary(codonTable)
    for AA in AAset:
        result.write('\"' + AA + '\",')
        for i in range(len(codonList[AA])):
            if not i == 0:
                result.write(',')
            result.write('\"' + codonList[AA][i] + '\",'+ str(codons[codonList[AA][i]]) + '\n')
    result.write('\"Total\",,' + str(totalAA))
    result.close()
    return



def _main():
    _console()
    input('Press any key to exit.')
    return



def _console():
    inputstr = input("Input functionality: ").strip().lower()
    if inputstr == 'help':
        print('possible operations:\n cleaning\n reverse\n rcomplement\n transcript\n rtranscript\n translate\n settings\n')
        _console()
        return

    elif inputstr == 'cleaning':
        print('Processing data...')
        data = fasta_reader('AA')
        fasta_writer(data,100,False)
        print('Task finished.\n')

    elif inputstr == 'reverse':
        print('Processing data...')
        data = fasta_reader('DNAX')
        processed = []
        for item in data:          
            processed.append((item[0], reverse(item[1])))
        fasta_writer(processed)
        print('Task finished.')

    elif inputstr == 'rcomplement':
        print('Processing data...')
        data = fasta_reader('DNA')
        processed = []
        for item in data:
            processed.append((item[0], reverse_complement(item[1])))
        fasta_writer(processed)
        print('Task finished.')

    elif inputstr == 'transcript':
        print('Processing data...')
        data = fasta_reader('DNA')
        processed = []
        for item in data:
            processed.append((item[0], transcript(item[1])))
        fasta_writer(processed)
        print('Task finished.')

    elif inputstr == 'rtranscript':
        print('Processing data...')
        data = fasta_reader('RNA')
        processed = []
        for item in data:
            processed.append((item[0], reverse_transcript(item[1])))
        fasta_writer(processed)
        print('Task finished.')    

    elif inputstr == 'translate':
        print('Processing data...')
        data = fasta_reader('DNA')
        processed = translate(data)
        fasta_writer(processed)
        print('Task finished.')
    else:
        print('Error: Uninterpretable input!')
        _console()
    return


if __name__ == '__main__':
    # execute only if run as a script
    _main()
