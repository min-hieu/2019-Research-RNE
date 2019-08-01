from strkernel.mismatch_kernel import MismatchKernel
from strkernel.mismatch_kernel import preprocess
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
import numpy as np
from numpy import random
import pandas as pd



# load the data
posSeq = [seq.seq for seq in SeqIO.parse('testing.fasta', 'fasta')]
negSeq = [seq.seq for seq in SeqIO.parse('testing.fasta', 'fasta')]

posX = preprocess(random.choice(posSeq, 2000))
negX = preprocess(random.choice(negSeq, 2000))

# label positive data as 1, negative as 0
posY = np.ones(len(posX), dtype=int)
negY = np.zeros(len(negX), dtype=int)

#compute mismatch kernels used in subsequent SVM training and testing
posKernels = MismatchKernel(l=4, k=5, m=1).get_kernel(posX).kernel
negKernels = MismatchKernel(l=4, k=5, m=1).get_kernel(negX).kernel

# merge data
X = np.concatenate([posKernels,negKernels])
y = np.concatenate([posY,negY])

# split training and test data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.30)

clf = SVC()
clf.fit(X_train, y_train)
