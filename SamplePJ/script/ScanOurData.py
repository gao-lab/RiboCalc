import pdb

import pandas as pd
import numpy as np
from sklearn import preprocessing
import scipy.stats as stats
import keras
np.random.seed(1337)
from pysam import Fastafile

import glob
from keras.preprocessing import sequence
from keras.optimizers import RMSprop
from keras.models import Sequential
from keras.layers.core import Dense
from keras.layers.core import Dropout
from keras.layers.core import Activation
from keras.layers.core import Flatten
from keras.layers.convolutional import Conv1D
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns

def one_hot_encode(seq, seq_len=50):
    # Dictionary returning one-hot encoding of nucleotides.
    nuc_d = {'a': [1, 0, 0, 0], 'c': [0, 1, 0, 0], 'g': [0, 0, 1, 0], 't': [0, 0, 0, 1], 'n': [0, 0, 0, 0]}

    # Creat empty matrix.
    vectors = np.empty([seq_len, 4])

    # Iterate through UTRs and one-hot encode
    seq = seq.lower()
    a = np.array([nuc_d[x] for x in seq])
    vectors = a
    return vectors


def GeneRate_model(border_mode='same', inp_len=50, nodes=40, layers=3, filter_len=8, nbr_filters=120,
                dropout1=0, dropout2=0.0, dropout3=0.0, nb_epoch=3):
    ''' Build model archicture and fit.'''
    model = Sequential()
    if layers >= 1:
        model.add(Conv1D(activation="relu", input_shape=(inp_len, 4), padding=border_mode, filters=nbr_filters,
                         kernel_size=filter_len))
    if layers >= 2:
        model.add(Conv1D(activation="relu", input_shape=(inp_len, 1), padding=border_mode, filters=nbr_filters,
                         kernel_size=filter_len))
        model.add(Dropout(dropout1))
    if layers >= 3:
        model.add(Conv1D(activation="relu", input_shape=(inp_len, 1), padding=border_mode, filters=nbr_filters,
                         kernel_size=filter_len))
        model.add(Dropout(dropout2))
    model.add(Flatten())

    model.add(Dense(nodes))
    model.add(Activation('relu'))
    model.add(Dropout(dropout3))

    model.add(Dense(1))
    model.add(Activation('linear'))

    # compile the model
    adam = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
    model.compile(loss='mean_squared_error', optimizer=adam)

    return model


def test_data(trainlabel, model, test_seq):
    '''Predict mean ribosome load using model and test set UTRs'''

    # Scale the test set mean ribosome load
    scaler = preprocessing.StandardScaler()
    scaler.fit(np.asarray(trainlabel).reshape(-1, 1))

    # Make predictions
    predictions = model.predict(test_seq).reshape(-1)

    return scaler.inverse_transform(predictions)
def r2(x,y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    return r_value**2


def PairData(Fasta, RiboTPM):
    """

    :param Fasta:
    :param RiboTPM:
    :return:
    """
    seqlist = []
    labellist = []
    orilabel = []

    f = open(RiboTPM,"r").readlines()
    for line in f:

        tem = line.replace("\n","").split("\t")
        seq = Fasta.fetch(tem[0])
        label = float(tem[1])
        seqtemlist, labeltemlist = GeneRateModelInput(seq, label, seqlen=50)
        orilabel.append(label)
        seqlist = seqlist + seqtemlist
        labellist = labellist + labeltemlist

    return np.asarray(seqlist), np.asarray(labellist), orilabel


def GeneRateModelInput(seq, label, seqlen=50):
    """

    :param seqlist:
    :param labellist:
    :param seqlen:
    :return:
    """
    SeqRealLen = len(seq)
    if SeqRealLen <= seqlen:
        seq = seq+"N"*(seqlen-SeqRealLen)
        SeqOutlist = [one_hot_encode(seq)]
        labelOutlist = [label]

    else:
        SeqOutlist = []
        labelOutlist = []
        for i in range(SeqRealLen-seqlen+1):
            sequse = seq[i:i+seqlen]
            SeqOutlist.append(one_hot_encode(sequse))
            labelOutlist.append(label)

    return SeqOutlist,labelOutlist



def loadData(path):
    """

    :param path:
    :return:
    """
    TrainSeqPath = path+"/coding_features.train.fasta"
    TrainRiboTPM = path+"/coding_features.train.RiboTPM"
    Allfasta = Fastafile(TrainSeqPath)

    TrainSeq, Trainlabel,OriTrainlabel = PairData(Allfasta, TrainRiboTPM)

    TestSeqPath = path+"/coding_features.test.fasta"
    TestRiboTPM = path+"/coding_features.test.RiboTPM"
    Allfasta = Fastafile(TestSeqPath)

    TestSeq, Testlabel,OriTestlabel = PairData(Allfasta, TestRiboTPM)

    return TrainSeq, Trainlabel, TestSeq, Testlabel,OriTrainlabel


def draw(predictions,Testlabel, outputpath,r):
    """

    Returns:

    """
    plt.scatter(predictions, Testlabel)
    plt.ylabel("Real RiboTPM")
    plt.xlabel("Model predict value")
    plt.text(6,1, "r-squared: "+str(round(r,4)))

    plt.savefig(outputpath+"/test.png")
    plt.close()


def mkdir(path):
    """
    Determine if the path exists, if it does not exist, generate this path
    :param path: Path to be generated
    :return:
    """
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return (True)
    else:
        return (False)

if __name__ == '__main__':
    import os
    os.environ["CUDA_VISIBLE_DEVICES"] = "6"

    datapath = "../raw_data/"

    TrainSeq, Trainlabel, TestSeq, Testlabel,OriTrainlabel = loadData(datapath)



    model = GeneRate_model(nb_epoch=3,border_mode='same',
                       inp_len=50, nodes=40, layers=3, nbr_filters=120, filter_len=8, dropout1=0,
                       dropout2=0, dropout3=0.2)

    model = keras.models.load_model('/rd1/user/lijy/human_5utr_modeling/modeling/saved_models/main_MRL_model.hdf5')


    predictions = test_data(OriTrainlabel, model=model, test_seq=TestSeq)
    r = r2(Testlabel, predictions)
    print('r-squared = ', r)
    outputpath = "../result/"
    mkdir(outputpath)
    np.savetxt("../result/predictions.txt",predictions)
    np.savetxt("../result/Testlabel.txt",Testlabel)
    outoutCSV = pd.DataFrame({"predictions":predictions, "Testlabel":Testlabel})
    outoutCSV.to_csv(outputpath+"/result.csv")
    draw(predictions,Testlabel, outputpath,r)






