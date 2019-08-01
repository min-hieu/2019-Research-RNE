# -*- coding: utf-8 -*-
# Python 3.7
# author Hyunjong Byun, Hieu
# Form implementation generated from reading ui file 'C:\Users\Admin\Desktop\PrimerBuilder\PrimerBuilder.ui'
#
# WARNING! All changes made in this file will be lost!
# Version 0.1()

import sys
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QLineEdit, QFileDialog, QInputDialog
from PyQt5.QtGui import QIcon, QGuiApplication
from PyQt5.Qt import QPushButton, QTextEdit
from SequenceProcessor import DNAset, DNAXset, RNAset, RNAXset
from SequenceProcessor import fasta_reader, fasta_writer
from SequenceProcessor import inverse_dictionary, mode_set
from SequenceProcessor import reverse_complement

enzset = {
    'XbaI': 'TCTAGA', 'NdeI': 'CATATG', 'KpnI': 'GGTACC', 'SacI': 'GAGCTC', \
    'HindIII': 'AAGCTT', 'EcoRI': 'GAATTC', 'SacII': 'CCGCGG'
}

tagset = {
    'His6': 'CATCACCATCACCATCAC', 'His6_alt': 'CACCACCACCACCACCAC'
}


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

        temp_seq = reverse_complement(temp_seq)
        if not enz == 'None':
            temp_seq = 'GGG' + reverse_complement(enzset[enz]) + temp_tag + temp_seq
    return (temp_header, temp_seq)


class Ui_PrimerBuilder(object):
    def setupUi(self, PrimerBuilder):
        self.onlyNum = QtGui.QIntValidator()
        # Layout and Spacers
        PrimerBuilder.setObjectName("PrimerBuilder")
        PrimerBuilder.resize(647, 580)
        PrimerBuilder.setMinimumSize(QtCore.QSize(647, 580))
        PrimerBuilder.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("H:/LOGO/ksa_logo_TWp_icon.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        PrimerBuilder.setWindowIcon(icon)
        self.centralwidget = QtWidgets.QWidget(PrimerBuilder)
        self.centralwidget.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setContentsMargins(-1, 11, -1, -1)
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")

        # label1
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setObjectName("label")
        self.horizontalLayout_10.addWidget(self.label)

        # radioDNA
        self.radio_DNA = QtWidgets.QRadioButton(self.centralwidget)
        self.radio_DNA.setObjectName("radio_DNA")
        self.buttonGroup_2 = QtWidgets.QButtonGroup(PrimerBuilder)
        self.buttonGroup_2.setObjectName("buttonGroup_2")
        self.buttonGroup_2.addButton(self.radio_DNA)
        self.horizontalLayout_10.addWidget(self.radio_DNA)
        self.radio_DNA.setChecked(True)

        # radioDNAX
        self.radio_DNAX = QtWidgets.QRadioButton(self.centralwidget)
        self.radio_DNAX.setObjectName("radio_DNAX")
        self.buttonGroup_2.addButton(self.radio_DNAX)
        self.horizontalLayout_10.addWidget(self.radio_DNAX)

        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout_10)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")

        # labelDNA
        self.label_DNA = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("MS Shell Dlg 2")
        self.label_DNA.setFont(font)
        self.label_DNA.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.label_DNA.setObjectName("label_DNA")
        self.horizontalLayout_2.addWidget(self.label_DNA)
        self.text_DNA = QtWidgets.QTextEdit(self.centralwidget)
        self.text_DNA.viewport().setProperty("cursor", QtGui.QCursor(QtCore.Qt.IBeamCursor))
        self.text_DNA.setPlaceholderText("Insert DNA sequence here!")
        self.text_DNA.setObjectName("text_DNA")

        # textDNA
        self.horizontalLayout_2.addWidget(self.text_DNA)

        # btnFileDNA
        self.btn_file_DNA = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_file_DNA.sizePolicy().hasHeightForWidth())
        self.btn_file_DNA.setSizePolicy(sizePolicy)
        self.btn_file_DNA.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.btn_file_DNA.setText("")
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("folder.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.btn_file_DNA.setIcon(icon1)
        self.btn_file_DNA.setObjectName("btn_file_DNA")
        self.horizontalLayout_2.addWidget(self.btn_file_DNA)
        self.btn_file_DNA.setShortcut("Ctrl+O")
        self.btn_file_DNA.clicked.connect(lambda: self.file_open())

        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setContentsMargins(0, 0, -1, -1)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.tm = QtWidgets.QLabel(self.centralwidget)
        self.tm.setAlignment(QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.tm.setObjectName("tm1")
        self.horizontalLayout_3.addWidget(self.tm)
        spacerItem1 = QtWidgets.QSpacerItem(70, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem1)

        # textTM
        self.text_tm = QtWidgets.QLineEdit(self.centralwidget)
        self.text_tm.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.text_tm.sizePolicy().hasHeightForWidth())
        self.text_tm.setSizePolicy(sizePolicy)
        self.text_tm.setMinimumSize(QtCore.QSize(12, 26))
        self.text_tm.setMaximumSize(QtCore.QSize(30, 27))
        self.text_tm.setWhatsThis("")
        self.text_tm.setObjectName("text_tm")
        self.horizontalLayout_3.addWidget(self.text_tm)
        self.text_tm.setValidator(self.onlyNum)
        self.text_tm.returnPressed.connect(lambda: self._Tm_select(self.text_tm.text()))

        # tm_wrn
        self.tm_wrn = QtWidgets.QLabel(self.centralwidget)
        self.tm_wrn.setObjectName("tm_wrn")
        self.horizontalLayout_3.addWidget(self.tm_wrn)
        self.tm_wrn.setStyleSheet('color: red')

        spacerItem2 = QtWidgets.QSpacerItem(285, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem2)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setContentsMargins(-1, 1, -1, -1)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.ResEnzyme1 = QtWidgets.QLabel(self.centralwidget)
        self.ResEnzyme1.setAlignment(QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.ResEnzyme1.setObjectName("ResEnzyme1")
        self.horizontalLayout_4.addWidget(self.ResEnzyme1)
        spacerItem3 = QtWidgets.QSpacerItem(9, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem3)

        # text_re1
        self.text_re1 = QtWidgets.QTextEdit(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.text_re1.sizePolicy().hasHeightForWidth())
        self.text_re1.setSizePolicy(sizePolicy)
        self.text_re1.setMinimumSize(QtCore.QSize(235, 26))
        self.text_re1.setMaximumSize(QtCore.QSize(67, 26))
        self.text_re1.viewport().setProperty("cursor", QtGui.QCursor(QtCore.Qt.IBeamCursor))
        self.text_re1.setObjectName("text_re1")
        self.horizontalLayout_4.addWidget(self.text_re1)

        self.re1_wrn = QtWidgets.QLabel(self.centralwidget)
        self.re1_wrn.setText("")
        self.re1_wrn.setObjectName("re1_wrn")
        self.horizontalLayout_4.addWidget(self.re1_wrn)

        spacerItem4 = QtWidgets.QSpacerItem(26, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem4)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.ResEnzyme2 = QtWidgets.QLabel(self.centralwidget)
        self.ResEnzyme2.setAlignment(QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.ResEnzyme2.setObjectName("ResEnzyme2")
        self.horizontalLayout_5.addWidget(self.ResEnzyme2)
        spacerItem5 = QtWidgets.QSpacerItem(9, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_5.addItem(spacerItem5)

        # text_re2
        self.text_re2 = QtWidgets.QTextEdit(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.text_re2.sizePolicy().hasHeightForWidth())
        self.text_re2.setSizePolicy(sizePolicy)
        self.text_re2.setMinimumSize(QtCore.QSize(67, 26))
        self.text_re2.setMaximumSize(QtCore.QSize(235, 28))
        self.text_re2.viewport().setProperty("cursor", QtGui.QCursor(QtCore.Qt.IBeamCursor))
        self.text_re2.setObjectName("text_re2")
        self.horizontalLayout_5.addWidget(self.text_re2)

        self.re2_wrn = QtWidgets.QLabel(self.centralwidget)
        self.re2_wrn.setText("")
        self.re2_wrn.setObjectName("re2_wrn")
        self.horizontalLayout_5.addWidget(self.re2_wrn)

        spacerItem6 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_5.addItem(spacerItem6)
        self.verticalLayout.addLayout(self.horizontalLayout_5)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        spacerItem7 = QtWidgets.QSpacerItem(26, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_7.addItem(spacerItem7)
        self.tag = QtWidgets.QLabel(self.centralwidget)
        self.tag.setAlignment(QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.tag.setObjectName("tag")
        self.horizontalLayout_7.addWidget(self.tag)

        self.wrn_tag = QtWidgets.QLabel(self.centralwidget)
        self.wrn_tag.setText("")
        self.wrn_tag.setObjectName("wrn_tag")
        self.horizontalLayout_7.addWidget(self.wrn_tag)
        self.tm_wrn.setStyleSheet('color: red')

        spacerItem8 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_7.addItem(spacerItem8)
        self.verticalLayout.addLayout(self.horizontalLayout_7)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")

        # tagNO
        self.tag_no = QtWidgets.QRadioButton(self.centralwidget)
        self.tag_no.setObjectName("tag_no")
        self.buttonGroup = QtWidgets.QButtonGroup(PrimerBuilder)
        self.buttonGroup.setObjectName("buttonGroup")
        self.buttonGroup.addButton(self.tag_no)
        self.horizontalLayout_6.addWidget(self.tag_no)
        self.tag_no.setChecked(True)

        spacerItem8 = QtWidgets.QSpacerItem(32, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_6.addItem(spacerItem8)

        # tagHIS
        self.tag_his = QtWidgets.QRadioButton(self.centralwidget)
        self.tag_his.setObjectName("tag_his")
        self.buttonGroup.addButton(self.tag_his)
        self.horizontalLayout_6.addWidget(self.tag_his)
        spacerItem9 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_6.addItem(spacerItem9)

        # tagOTHER
        self.tag_other = QtWidgets.QRadioButton(self.centralwidget)
        self.tag_other.setObjectName("tag_other")

        self.buttonGroup.addButton(self.tag_other)

        self.horizontalLayout_6.addWidget(self.tag_other)
        self.text_tag = QtWidgets.QTextEdit(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.text_tag.sizePolicy().hasHeightForWidth())
        self.text_tag.setSizePolicy(sizePolicy)
        self.text_tag.setMinimumSize(QtCore.QSize(0, 0))
        self.text_tag.setMaximumSize(QtCore.QSize(16777215, 26))
        self.text_tag.viewport().setProperty("cursor", QtGui.QCursor(QtCore.Qt.IBeamCursor))
        self.text_tag.setObjectName("text_tag")
        self.horizontalLayout_6.addWidget(self.text_tag)

        self.verticalLayout.addLayout(self.horizontalLayout_6)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout.setObjectName("horizontalLayout")

        self.sense = QtWidgets.QLabel(self.centralwidget)
        self.sense.setObjectName("sense")
        self.horizontalLayout.addWidget(self.sense)

        # senseYES
        self.sense_yes = QtWidgets.QCheckBox(self.centralwidget)
        self.sense_yes.setObjectName("sense_yes")
        self.horizontalLayout.addWidget(self.sense_yes)
        spacerItem10 = QtWidgets.QSpacerItem(24, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem10)

        # senseNO
        self.sense_anti = QtWidgets.QCheckBox(self.centralwidget)
        self.sense_anti.setObjectName("sense_anti")
        self.horizontalLayout.addWidget(self.sense_anti)
        spacerItem11 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem11)

        self.verticalLayout.addLayout(self.horizontalLayout)
        spacerItem12 = QtWidgets.QSpacerItem(18, 12, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.verticalLayout.addItem(spacerItem12)
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_8.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        spacerItem13 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem13)

        # BuildPrime
        self.buildPrime = QtWidgets.QPushButton(self.centralwidget)
        self.buildPrime.setMinimumSize(QtCore.QSize(110, 0))
        font = QtGui.QFont()
        font.setFamily("Comic Sans MS")
        font.setPointSize(9)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.buildPrime.setFont(font)
        self.buildPrime.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.buildPrime.setObjectName("buildPrime")
        self.horizontalLayout_8.addWidget(self.buildPrime)

        spacerItem14 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem14)
        self.verticalLayout.addLayout(self.horizontalLayout_8)
        spacerItem15 = QtWidgets.QSpacerItem(20, 12, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.verticalLayout.addItem(spacerItem15)

        self.result = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("MS Shell Dlg 2")
        font.setPointSize(10)
        self.result.setFont(font)
        self.result.setObjectName("result")
        self.verticalLayout.addWidget(self.result)

        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")

        # textRESULT
        self.result_text = QtWidgets.QTextEdit(self.centralwidget)
        self.result_text.setMinimumSize(QtCore.QSize(0, 15))
        self.result_text.viewport().setProperty("cursor", QtGui.QCursor(QtCore.Qt.IBeamCursor))
        self.result_text.setObjectName("result_text")
        self.horizontalLayout_9.addWidget(self.result_text)
        self.result_text.setStyleSheet('color: black')

        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setContentsMargins(0, -1, -1, -1)
        self.verticalLayout_2.setObjectName("verticalLayout_2")

        self.actionCopy = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.actionCopy.sizePolicy().hasHeightForWidth())
        self.actionCopy.setSizePolicy(sizePolicy)
        self.actionCopy.setMinimumSize(QtCore.QSize(0, 30))
        self.actionCopy.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.actionCopy.setText("")
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("copy.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionCopy.setIcon(icon2)
        self.actionCopy.setObjectName("actionCopy")
        self.verticalLayout_2.addWidget(self.actionCopy)

        self.actionErase = QtWidgets.QPushButton(self.centralwidget)
        self.actionErase.setMinimumSize(QtCore.QSize(0, 30))
        self.actionErase.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.actionErase.setText("")
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap("broom.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionErase.setIcon(icon3)
        self.actionErase.setObjectName("actionErase")
        self.verticalLayout_2.addWidget(self.actionErase)

        self.horizontalLayout_9.addLayout(self.verticalLayout_2)
        self.verticalLayout.addLayout(self.horizontalLayout_9)
        spacerItem16 = QtWidgets.QSpacerItem(13, 5, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.verticalLayout.addItem(spacerItem16)
        self.credit_1 = QtWidgets.QLabel(self.centralwidget)
        self.credit_1.setObjectName("credit_1")
        self.verticalLayout.addWidget(self.credit_1)
        self.credit_2 = QtWidgets.QLabel(self.centralwidget)
        self.credit_2.setObjectName("credit_2")
        self.verticalLayout.addWidget(self.credit_2)
        PrimerBuilder.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(PrimerBuilder)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 647, 26))
        self.menubar.setObjectName("menubar")
        PrimerBuilder.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(PrimerBuilder)
        self.statusbar.setObjectName("statusbar")
        PrimerBuilder.setStatusBar(self.statusbar)

        self.retranslateUi(PrimerBuilder)
        QtCore.QMetaObject.connectSlotsByName(PrimerBuilder)

        self.buildPrime.clicked.connect(lambda: self.buildNOW())

    @pyqtSlot()
    def _Tm_select(self, tm__):
        if tm__ == '':
            pass
        elif int(tm__) < 50:
            self.tm_wrn.setText('Warning: Too low Tm.')
        else:
            self.tm_wrn.setText('')

        return tm__

    def _primer_select(self, fr, primer__):
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

    def _tag_select(self, mode):
        tag = self.buttonGroup.checkedButton().text()[2:].strip()

        if tag == 'Other sequence:':
            tag = self.text_tag.toPlainText()
            comset = mode_set(mode)
            temp_tag = ''
            check_err = 0
            for char in tag:
                if char in comset:
                    temp_tag += char
                else:
                    check_err += 1

            err_tag = '    ' + str(check_err) + ' INVALID CHARACTER(s) WERE IGNORED'
            tag = temp_tag
            self.wrn_tag.setText(err_tag)

        return tag

    def _dna_select(self, dna, mode):

        lines = dna.split('\n')
        comset = mode_set(mode)
        out = []
        started = False
        temp_header = ''
        temp_seq = []
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
        out.append((temp_header, temp_seq))

        print(lines)

        return out

    def buildNOW(self):
        mode = self.buttonGroup_2.checkedButton().text()
        forward = self._primer_select('f', self.text_re1.toPlainText())
        reverse = self._primer_select('r', self.text_re2.toPlainText())
        r_data = self._dna_select(self.text_DNA.toPlainText(), mode)
        r_tm = self._Tm_select(self.text_tm.text())
        r_tag = self._tag_select(mode)

        for i in [forward, reverse, r_data, r_tm, r_tag]:
            if i == 'None':
                self.result_err_text += 'Please fill in the information. \n'

                if self.tag_no.isChecked():
                    if self.sense_yes.isChecked() == False and self.sense_anti.isChecked() == False:
                        self.result_err_text += 'Please check sense or antisense'
                if i == 'err':
                    self.result_err_text += 'Error: Forward or backward primer does not exist in dictionary "enzset". \n'

                self.result_text.setText(self.result_err_text)

            else:
                if self.sense_yes.isChecked():
                    ftag = r_tag
                if self.sense_anti.isChecked():
                    rtag = r_tag
                processed = []
                for item in r_data:
                    processed.append(primer_builder('f', item, forward, ftag, r_tm, mode))
                    processed.append(primer_builder('r', item, reverse, rtag, r_tm, mode))
                    print(processed)
                fasta_writer(processed, 999999999, False)
                print('Task finished.')
                f_ = open("./output.txt", 'r')

                with f_:
                    text = f_.read()
                    self.result_text.setText(text)


    def result_err(self, err__):

        self.result_text.setStyleSheet('color: red')
        self.result_text.setText(err__)

    def file_open(self):
        filename = QFileDialog.getOpenFileName(None, 'Select Fasta File', '.txt')
        file = open(list(filename)[0], 'r')

        with file:
            text = file.read()
            self.text_DNA.setText(text)

    def retranslateUi(self, PrimerBuilder):
        _translate = QtCore.QCoreApplication.translate
        PrimerBuilder.setWindowTitle(_translate("PrimerBuilder", "PrimerBuilder"))
        self.label.setText(_translate("PrimerBuilder", "Input type:"))
        self.radio_DNA.setText(_translate("PrimerBuilder", "DNA"))
        self.radio_DNAX.setText(_translate("PrimerBuilder", "DNAX"))
        self.label_DNA.setText(_translate("PrimerBuilder", "DNA SEQUENCE"))
        self.tm.setText(_translate("PrimerBuilder", "Tm"))
        self.ResEnzyme1.setText(_translate("PrimerBuilder", "Res Enzyme 1"))
        self.ResEnzyme2.setText(_translate("PrimerBuilder", "Res Enzyme 2"))
        self.tag.setText(_translate("PrimerBuilder", "Tag?"))
        self.tag_no.setText(_translate("PrimerBuilder", "No Tag"))
        self.tag_his.setText(_translate("PrimerBuilder", "His-Tag"))
        self.tag_other.setText(_translate("PrimerBuilder", "Other sequence:"))
        self.text_tag.setPlaceholderText(_translate("PrimerBuilder", "Insert custom tag here!"))
        self.sense.setText(_translate("PrimerBuilder", "Tag at sense and/or anti-sense? "))
        self.sense_yes.setText(_translate("PrimerBuilder", "sense"))
        self.sense_anti.setText(_translate("PrimerBuilder", "anti-sense"))
        self.buildPrime.setText(_translate("PrimerBuilder", "Build Primer"))
        self.result.setText(_translate("PrimerBuilder", "Result: "))
        self.credit_1.setText(_translate("PrimerBuilder",
                                         "PrimerBuilder Reborn v1.2 powered by SequenceProcessor v2.2 by Hyunjong Byun."))
        self.credit_2.setText(_translate("PrimerBuilder", "implemented by Minh Hieu (18-206)"))


def main():
    app = QtWidgets.QApplication(sys.argv)
    PrimerBuilder = QtWidgets.QMainWindow()
    ui = Ui_PrimerBuilder()
    ui.setupUi(PrimerBuilder)
    PrimerBuilder.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
