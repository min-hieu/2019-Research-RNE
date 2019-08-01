import sys
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
from PyQt5.uic import loadUi
from SequenceProcessor import mode_set, reverse_complement

enzset = {
    'XbaI': 'TCTAGA', 'NdeI': 'CATATG', 'KpnI': 'GGTACC', 'SacI': 'GAGCTC', \
    'HindIII': 'AAGCTT', 'EcoRI': 'GAATTC', 'SacII': 'CCGCGG'
}

tagset = {
    'His6': 'CATCACCATCACCATCAC', 'His6_alt': 'CACCACCACCACCACCAC'
}

# primer_builder mod by Hieu
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
                if calc_Tm >= int(Tm):
                    break
            if not enz == 'None':
                temp_seq = 'GGG' + enzset[enz] + temp_tag + temp_seq
        else:
            for i in range(len(item[1])):
                index = len(item[1]) - 1 - i
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

def _check_any(list, key):
    for i in list:
        if i == key:
            return True
    return False


class PrimerBuilder(QMainWindow):
    def __init__(self):
        super(PrimerBuilder, self).__init__()
        loadUi("PrimerBuilder.ui", self)
        self.button_file.clicked.connect(lambda: self.file_open())
        self.button_build.clicked.connect(lambda: self.buildNOW())
        self.onlyNum = QtGui.QIntValidator()
        self.text_tm.setValidator(self.onlyNum)

    @pyqtSlot()

    def _dna_select(dna, mode):
        lines = dna.split('\n')
        comset = mode_set(mode)
        out = []
        started = False
        temp_header = ''
        temp_seq_list = []

        for line in lines:
            line_lstrip = line.lstrip()
            if line_lstrip.strip() == '':
                continue
            elif line_lstrip[0] == '>':
                if started:
                    if temp_seq_list == []:
                        break
                    else:
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
                        print('Error: character \"' + l[i] + '\" does not match with specified type \"' \
                              + mode + '\" in function fasta_reader(' + dna + ').')
        temp_seq = ''.join(temp_seq_list)
        out.append([temp_header, temp_seq])

        print(lines)

        return out

    def _Tm_select(self, tm__):
        if int(tm__) < 50:
            self.popup_err('Warning: Too low Tm.')

        return tm__

    def _primer_select(self, fr, primer__):
        try:
            if fr == 'f':
                _text = 'Forward'
            elif fr == 'r':
                _text = 'Reverse'
            primer = primer__.strip()
            if primer == '' or primer == 'None':
                return 'None'
            elif primer not in enzset:
                print('Error: ' + _text + ' primer does not exist in dictionary "enzset".')
                primer = 'err'
            return primer
        except:
            self.popup_err("Check Res Enzyme Boxes")

    def _tag_select(self, mode):
        tag = self.buttons_tag.checkedButton().text()

        if tag == 'Other:':
            tag = self.text_tag.toPlainText()
            comset = mode_set(mode)
            temp_tag = ''
            check_err = 0
            for char in tag:
                if char in comset:
                    temp_tag += char
                else:
                    check_err += 1

            self.popup_err('    ' + str(check_err) + ' INVALID CHARACTER(s) WERE IGNORED IN CUSTOM TAG')
            tag = temp_tag

        elif tag == 'His-Tag':
            tag = 'His6'
        else:
            tag = ''

        return tag

    def buildNOW(self):
        mode = self.buttons_dna.checkedButton().text()
        r_data = _dna_select(self.text_dna.toPlainText(), mode)
        forward = self._primer_select('f', self.text_res1.toPlainText())
        reverse = self._primer_select('r', self.text_res2.toPlainText())
        r_tm = self._Tm_select(self.text_tm.text())
        r_tag = self._tag_select(mode)

        try:
            result_err_text = ''
            if _check_any([forward, reverse, r_data, r_tm, r_tag], 'None'):
                result_err_text += 'Please fill in the information! \n'

                if self.radio_noTag.isChecked():
                    if self.check_sense.isChecked() == False and self.check_antiSense.isChecked() == False:
                        result_err_text += 'Please check sense or antisense \n'
                if _check_any([forward, reverse, r_data, r_tm, r_tag], 'err'):
                    self.result_err_text += 'Forward or backward primer does not exist in dictionary "enzset". \n'

                self.popup_err(result_err_text)

            else:

                ftag = rtag = 'None'
                if self.check_sense.isChecked():
                    ftag = r_tag
                if self.check_antiSense.isChecked():
                    rtag = r_tag
                processed = []
                for item in r_data:
                    processed.append(primer_builder('f', item, forward, ftag, r_tm, mode))
                    # processed.append(primer_builder('r', item, reverse, rtag, r_tm, mode))
                    print(processed)
                self.hieu_fasta_writer(processed)
        except:
            pass

    def popup_err(self, err__):
        msgBox = QtWidgets.QMessageBox()
        msgBox.setIcon(QtWidgets.QMessageBox.Warning)
        msgBox.setWindowTitle("Error")
        erricon = QtGui.QIcon()
        erricon.addPixmap(QtGui.QPixmap("error.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        msgBox.setWindowIcon(erricon)
        msgBox.setText(err__)
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        msgBox.exec_()

    def hieu_fasta_writer(self, list):
        with open('output.txt', 'w') as f:
            for i in list:
                f.write(i + '\n')

        with open('output.txt', 'r') as f:
            text = f.read()
            self.text_result.setText(text)

    def file_open(self):
        try:
            filename = QFileDialog.getOpenFileName(None, 'Select Fasta File', '.txt')
            file = open(list(filename)[0], 'r')

            with file:
                text = file.read()
                self.text_dna.setText(text)
        except:
            self.popup_err("wrong file format")


def main():
    app = QApplication(sys.argv)
    widget = PrimerBuilder()
    widget.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
