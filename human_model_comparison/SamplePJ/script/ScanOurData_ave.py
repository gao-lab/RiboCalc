import pdb

import pandas as pd
import numpy as np
from sklearn import preprocessing
import scipy.stats as stats
np.random.seed(1337)

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns



def r2(x,y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    # return r_value**2
    return r_value



def draw(predictions,Testlabel, outputpath,r):
    """

    Returns:

    """
    font2 = {'family': 'arial',
             'weight': 'normal',
             'size': 12,
             }
    plt.scatter(predictions, Testlabel,s=1)
    plt.ylabel("Observed Ribo-TPM",font2)
    plt.xlabel("Sample's model predicted Ribo-TPM",font2)
    # plt.text(4,1, "r-squared: "+str(round(r,4)))
    plt.text(-1,10, "r: "+str(round(r,4)))

    plt.savefig(outputpath)
    plt.savefig(outputpath.replace("png","eps"),format="eps")
    plt.close()

def PariAllData(predictionsALL,TestlabelALL):
    """

    :param predictionsALL:
    :param TestlabelALL:
    :return:
    """
    InitLabel = TestlabelALL[0]
    Testlabel = [InitLabel]
    predictions= []
    pretem = []
    for i in range(len(TestlabelALL)):
        if TestlabelALL[i]==InitLabel:
            pretem.append(predictionsALL[i])
        else:
            InitLabel = TestlabelALL[i]
            Testlabel.append(InitLabel)
            if i !=0:
                predictions.append(np.mean(pretem))
                pretem = [predictionsALL[i]]

        if i == len(TestlabelALL)-1:

            predictions.append(np.mean(pretem))


    return Testlabel, predictions



if __name__ == '__main__':
    outputpath = "../result/ave.png"
    predictionsALL = np.loadtxt("../result/predictions.txt",)
    TestlabelALL = np.loadtxt("../result/Testlabel.txt")

    Testlabel, predictions = PariAllData(predictionsALL, TestlabelALL)
    r = r2(Testlabel, predictions)
    outoutCSV = pd.DataFrame({"predictions":predictions, "Testlabel":Testlabel})
    outoutCSV.to_csv("../result/result_ave.csv")
    draw(predictions,Testlabel, outputpath,r)






