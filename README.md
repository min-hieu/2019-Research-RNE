# 2019-Research-R-E

This repo contain all of Hieu's script and Supercharging Algorithm.

### Inside **Main.py**, you will find:

#### •  *Super_seq* class
  This class is made to classify all the sequence as an object with features (sequence, charge, file location,etc) and functions related to supercharging the protein.
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Parameters

##### &nbsp;&nbsp;&nbsp; *seq* a str format of the sequence of the protein
##### &nbsp;&nbsp;&nbsp; *__raw_seq* raw format of the sequence (before formatting)
##### &nbsp;&nbsp;&nbsp; *charge* a list that contain the charge of every 20 amino acid
##### &nbsp;&nbsp;&nbsp; *consurf* (work in-progress) read the consurf scoring file obtained from the server of the protein
##### &nbsp;&nbsp;&nbsp; *name* name of the protein ex) TNFb
##### &nbsp;&nbsp;&nbsp; *location* location of the file containing the protein (default "")

#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Functions
##### &nbsp;&nbsp;&nbsp; *format_seq* format the input string of sequence
##### &nbsp;&nbsp;&nbsp; *seq_charge* obtain a list of the total charge of every 20 amino acid from the input sequence
##### &nbsp;&nbsp;&nbsp; *consurf_reader* read the scoring of the input .grade file of the protein
##### &nbsp;&nbsp;&nbsp; *display_lcd* draw the linear charge density graph of the protein
##### &nbsp;&nbsp;&nbsp; *compare_lcd* draw the linear charge density graph of the protein and another protein

#### &nbsp;&nbsp;&nbsp; ---------------------------------------------

#### •  Functions outside *Super_seq* class
##### &nbsp;&nbsp;&nbsp; *binding_site_converter* convert input binding site into a list ex) "10-14+17+19+21-23" -> [10,11,12,13,14,17,19,21,22,23]
##### &nbsp;&nbsp;&nbsp; the rest of the function is still in early development stages
