#!/usr/bin/python
# -*- coding: UTF-8 -*-
from rklib import utils,compress
import sys
#import log
import numpy as np
import scipy as sp
from scipy import stats
from rblib.mutilstats import SampleInfo,MatrixAnno,centring,normalize,c_permutation
import os
from rblib import gwaspls,statplot
class plsregress_opt(object):
    def __init__(self):
        self.log2tr = 0 
        self.phenotype = ""
        self.outdir = ""
        self.sampleinfo = ""
        self.prefix = ""
        self.fc = 0.000
        self.class_info = ""
        self.regmd = "pearsonr"
        self.rho = 0.000
        self.nlvs = 10
        self.do_predict = 0
        self.nfold = 4
        self.viptest = 0
        self.nviptest = 1000
        self.vipp = 0.05
        self.mcuve = 0
        self.mc_nlvs = 5
        self.mcratio = 0.5
        self.mctimes = 200
        self.nlvsoptimize = 0
        self.nlvsmax = 20
        self.mcvip = 0
        self.vip_mcvip = 0
        self.ndimentions = 2
        self.scoreplot = ""
        self.plot = 0
        self.matrix_anno = ""
        self.ex_optimize = 0
        self.ex_predict  = 0
        self.ex_sampleinfo = None
        self.ex_matrix = None
        #self.logger=log.initlog(os.getcwd()+os.path.sep+"plsregress.log")
    def destroy(self):
        #log.closelog()
        return 0

def __get_fc(class_array,data,traits,p):
    fc = []
    cidx_1 = np.where(traits == class_array[0])
    cidx_2 = np.where(traits == class_array[1])
    if (class_array[0] in traits) and (class_array[1] in traits):pass
    else:
        return None
    for i in range(p):
        tmpfc = np.mean(data[cidx_1,i]) - np.mean(data[cidx_2,i])
        fc.append(tmpfc)
    return fc

def __get_rho(regmd,data,traits,p):
    if regmd == "pearsonr":
        __regr_core = stats.pearsonr
    elif regmd == "spearmanr":
        __regr_core = stats.spearmanr
    else:
        return None
    rho = []
    for i in range(p):
        correlation,pvalue = __regr_core(data[:,i],traits)
        try:
            rho.append(correlation[0])
        except:
            rho.append(correlation)
    return rho

def check_dimension(p,nlvs,opt):
    if p < nlvs:
        sys.stderr.write("number of variables (p) < nLVs")
        return 1
    return 0

def viptest_output(fc_rho,pvalue,opt):
    ret_pvalue = None
    if opt.phenotype == "quantificat":
        cutoff = opt.rho
        xlabel = "Correlation coefficient"
    else:
        cutoff = opt.fc
        xlabel = "Fold change"
    fselect_stat = open(opt.outdir+os.path.sep+opt.prefix+".viptest.tab","w")
    ftatol_stat = open(opt.outdir+os.path.sep+opt.prefix+".viptest.total.tab","w")
    prefix,suffix = utils.parse_filename(os.path.basename(opt.matrix_anno))
    fselect_data = open(opt.outdir+os.path.sep+prefix+".viptest.anno","w")
    data = MatrixAnno()
    ret = data.parse_matrix_anno(opt.matrix_anno)
    tmp = []
    for i in range(data.p):
        anno = data.anno[i]
        tmpfcrho = fc_rho[i]
        tmppvalue = pvalue[i]
        tmpdata = np.asarray(data.data[:,i].T)[0].tolist()
        strdata = [format(i,".17g") for i in tmpdata]
        ftatol_stat.write("%s\t%.7g\t%.7g\n"%(anno,tmpfcrho,tmppvalue))
        if np.abs(tmpfcrho) > cutoff and tmppvalue < opt.vipp:
            tmp.append(tmppvalue)
            fselect_stat.write("%s\t%.7g\t%.7g\n"%(anno,tmpfcrho,tmppvalue))
            fselect_data.write(anno+"\t"+"\t".join(strdata)+"\n")
    fselect_data.close()
    fselect_stat.close()
    ftatol_stat.close()
    fig_prefix = opt.outdir+os.path.sep+opt.prefix+".viptest_vacoplot"
    #print pvalue
    statplot.vaco_plot(np.asarray(fc_rho),np.asarray(-1*np.log10(pvalue)),cutoff,-1*np.log10(opt.vipp),fig_prefix,xlabel,"$-log_{10}(probability)$",title=None)
    if tmp:
        ret_pvalue = tmp
    return ret_pvalue

def mcmd_output(fc_rho,stability,opt,best_stat):
    if opt.mcuve: ftype = ".viptest.mcuve."
    elif opt.mcvip:ftype = ".viptest.mcvip."
    #fselect_stat = open(opt.outdir+os.path.sep+opt.prefix+ftype+"tab","w")
    ftatol_stat = open(opt.outdir+os.path.sep+opt.prefix+ftype+"total.tab","w")
    ftatol_stat.write("#Best stat:%.5f\n"%best_stat)
    #prefix,suffix = utils.parse_filename(os.path.basename(opt.matrix_anno))
    #fselect_data = open(opt.outdir+os.path.sep+prefix+ftype+"anno","w")
    data = MatrixAnno()
    ret = data.parse_matrix_anno(opt.matrix_anno)
    for i in range(data.p):
        anno = data.anno[i]
        tmpfcrho = fc_rho[i]
        #tmppvalue = ret_pvalue[i]
        tmpstat = stability[i]
        tmpdata = np.asarray(data.data[:,i].T)[0].tolist()
        strdata = [format(i,".17g") for i in tmpdata]
        #ftatol_stat.write("%s\t%.7g\t%.7g\t%.7g\n"%(anno,tmpfcrho,tmppvalue,tmpstat))
        ftatol_stat.write("%s\t%.7g\t%.7g\n"%(anno,tmpfcrho,tmpstat))
        #if tmpstat >=best_stat:
            #fselect_stat.write("%s\t%.7g\t%.7g\t%.7g\n"%(anno,tmpfcrho,tmppvalue,tmpstat))
            #fselect_data.write(anno+"\t"+"\t".join(strdata)+"\n")
    #fselect_stat.close()
    ftatol_stat.close()
    #fselect_data.close()
    return 0

def calVIPdat(X,Y,opt,p):
    if check_dimension(p,opt.nlvs,opt):return None
    result = gwaspls.plsgwas(X,Y,opt.nlvs)
    VIPmat = gwaspls.getVIP(result,Y)
    print("####VIPlist#####")
    print(VIPmat)
    return 0

def do_viptest(X,Y,opt,p):
    if check_dimension(p,opt.nlvs,opt):return None
    result = gwaspls.plsgwas(X,Y,opt.nlvs)
    VIPmat = gwaspls.getVIP(result,Y)
    print("####VIPlist#####")
    print(VIPmat)
    fname_vip = os.path.sep.join([opt.outdir,opt.prefix+".vipmat.dat"])
    fh_vip = open(fname_vip,'wb')
    h = gwaspls.VIPoutput(VIPmat,fh_vip,p)
    if h == 0:sys.stderr.write("VIP Statics output to file %s"%fname_vip)
    else:
        sys.stderr.write("VIP Statics output failed")
        return None
    fh_vip.close()
    fname_per = os.path.sep.join([opt.outdir,opt.prefix+".permat.dat"])
    fh_per = open(fname_per,'wb+')
    h = gwaspls.VIPtest(X,Y,fh_per,opt.nlvs,ntimes=opt.nviptest)
    if h == 0: sys.stderr.write("VIP permutation Statics output to file %s"%fname_per)
    else:
        sys.stderr.write("VIP permutation Statics output failed")
        return None
    fh_per.close()
    pvalue = c_permutation(fname_vip,fname_per,p,opt.nviptest)
    if pvalue == None:
        sys.stderr.write("permutation test failed")
        return None
    else:sys.stderr.write("permutation test success")
    pvalue = np.asarray(pvalue)
    return pvalue

def stability_optimize(X,Y,opt,stability):
    step = 0.1
    s_range = np.arange(np.min(stability),np.max(stability),0.1)
    s_range = s_range[::-1]
    num = len(s_range)
    sys.stderr.write("stability calculation progress:\n")
    risks = []
    numbers = []
    svector = []
    flag = 0
    for i in range(num):
        s_tmp = s_range[i]
        idx = stability >= s_tmp
        tmp_numvars = np.sum(idx)
        if tmp_numvars > opt.mc_nlvs:
            flag = 1
            numbers.append(tmp_numvars)
            svector.append(s_tmp)
            risk = [0.0]*30
            for j in range(30):
                tmprisk = gwaspls.plscv(X[:,idx],Y,opt.mc_nlvs,opt.nfold,method=opt.phenotype)
                risk[j] = tmprisk
            risks.append(risk)
        sys.stderr.write("..... %s%%\n"%(str(np.round(float(i+1)/num*100,1))))
            #print s_tmp,risk
    sys.stderr.write("stability calculation Done!\n")
    if not risks:
        if flag == 0:
            sys.stderr.write("number of variables (p) < nLVs")
        return None
    mean_risk = np.mean(np.asarray(risks),axis=1)
    #std_risk  = np.std(np.asarray(risks),ddof=1,axis=1)
    numbers = np.asarray(numbers)
    if opt.mcuve:
        fig_prefix = opt.outdir+os.path.sep+opt.prefix+".mcuve_stability_optimize"
    elif opt.mcvip:
        fig_prefix = opt.outdir+os.path.sep+opt.prefix+".mcvip_stability_optimize"
    if opt.phenotype == "category": 
        ylabel1 = "Accuracy of cross validation"
    else:
        ylabel1 = "Root mean squares error ofcross validation"
    ylabel2 = "Numbers of selected features"
    if opt.phenotype == "category":
        up_risk = np.max(np.asarray(risks),axis=1)
        low_risk = np.min(np.asarray(risks),axis=1)
        best_risk = np.max(mean_risk)
    else:
        std_risk  = np.std(np.asarray(risks),ddof=1,axis=1)
        up_risk = mean_risk + std_risk
        low_risk = mean_risk - std_risk
        best_risk = np.min(mean_risk)
    risk_plot =  np.asarray((mean_risk,up_risk,low_risk))
    idx = mean_risk.tolist().index(best_risk)
    best_stat = svector[idx]
    title = "Best stability cutoff is %.6g"%best_stat
    ret = statplot.plotyy(np.asarray(svector),risk_plot,numbers,fig_prefix,"Stability",ylabel1,ylabel2,title)
    if ret == 0:
        return best_stat
    else:
        return None

def plsregress_script(opt,matrix_anno):
    ret = 1
    sampleinfo = SampleInfo()
    ret = sampleinfo.parse_sampleinfo(opt.sampleinfo,opt.phenotype)
    data = MatrixAnno()
    ret = data.parse_matrix_anno(matrix_anno)
    #print data.data
    if opt.phenotype == "quantificat":
        fc_rho = __get_rho(opt.regmd,data.data,sampleinfo.traits,data.p)
    elif opt.phenotype == "category":
        flag = 0
        class_infos = opt.class_info.split("v")
        if len(class_infos) !=2: flag=1
        try: class_infos = list(map(int,class_infos))
        except: flag=1
        if flag:
            sys.stderr.write("category compare failed, please check '-c %s'"%opt.class_info)
            return flag
        fc_rho = __get_fc(class_infos,data.data,sampleinfo.traits,data.p)
        if None == fc_rho:
            sys.stderr.write("category compare failed, please check '-c %s'"%opt.class_info)
            return 1
    #print fc_rho
    trainmean = centring(data.data)
    trainstd = normalize(data.data)
    ret_traits_mean = centring(sampleinfo.traits)
    ret_traits_std = 1
    if opt.phenotype != "category":
        ret_traits_std = normalize(sampleinfo.traits)
    else:
        sampleinfo.traits = np.sign(sampleinfo.traits)
    ## to do viptest
    flag = opt.viptest + opt.vip_mcvip + opt.mcuve + opt.mcvip
    if flag > 1:
        sys.stderr.write("viptest, mcuve, mcvip, vip_mcvip should choose one of them")
        return 1
    if opt.do_predict:## predict can shu 
        risk_cal,ypredict_cal = gwaspls.pls_calibration_risk(data.data,sampleinfo.traits,opt.nlvs,method=opt.phenotype,retpred=1)
        ypredict_cal = ypredict_cal * ret_traits_std + ret_traits_mean
        fout_predict = open("ypredict_cal.xls","w")
        for i in range(ypredict_cal.shape[0]):
            fout_predict.write(str(ypredict_cal[i,0])+"\n")
        fout_predict.close()
    if opt.ex_predict or opt.ex_optimize:
        exsampleinfo = SampleInfo()
        exsampleinfo.parse_sampleinfo(opt.ex_sampleinfo,opt.phenotype)
        exdata = MatrixAnno()
        exdata.parse_matrix_anno(opt.ex_matrix)
        ## do to get the model 
        if opt.ex_predict:
            ex_risk,ex_predict = gwaspls.pls_validataion_risk(data.data,sampleinfo.traits,trainmean,trainstd,ret_traits_mean,ret_traits_std,exdata.data,exsampleinfo.traits,opt.nlvs,method=opt.phenotype,retpred=1)
            fout_predict = open("ypredict_external.xls","w")
            fout_predict.write("#risk is : %.3f\n"%ex_risk)
            for i in range(ex_predict.shape[0]):
                fout_predict.write(str(ex_predict[i,0])+"\n")
            fout_predict.close()
        if opt.ex_optimize:
            risks = np.zeros((opt.nlvsmax,30))
            cal_risks = np.zeros(opt.nlvsmax)
            ex_risks = np.zeros(opt.nlvsmax)
            for i in range(1,opt.nlvsmax+1):
                for j in range(30):
                    risk = gwaspls.plscv(data.data,sampleinfo.traits,i,opt.nfold,method=opt.phenotype)
                    risks[i-1,j] = risk
                risk_cal = gwaspls.pls_calibration_risk(data.data,sampleinfo.traits,i,method=opt.phenotype)
                cal_risks[i-1] = risk_cal
                risks_ex = gwaspls.pls_validataion_risk(data.data,sampleinfo.traits,trainmean,trainstd,ret_traits_mean,ret_traits_std,exdata.data,exsampleinfo.traits,i,method=opt.phenotype)
                ex_risks[i-1] = risks_ex
                sys.stderr.write("nLVs: %d\n"%i)
            xlabel = "Number of latent variables"
            fig_prefix = opt.outdir+os.path.sep+opt.prefix+".nlvsoptimize.with_external"
            xticks_labels = list(map(str,range(1,opt.nlvsmax+1)))
            mean_risk = np.mean(np.asarray(risks),axis=1)
            if opt.phenotype == "category":
                ylabel = "Classification accuracy"
                up_risk = np.max(np.asarray(risks),axis=1)
                low_risk = np.min(np.asarray(risks),axis=1)
                risk_plot =  np.asarray((cal_risks,ex_risks,mean_risk,up_risk,low_risk))
                ret = statplot.plotline(np.arange(1,opt.nlvsmax+1),risk_plot*100,fig_prefix,xlabel,ylabel,['b-','g-','r-','r--','r--'],["Calibration","External Validation","Cross Validation"],title="Latentvariables selection",linewidth=1.5,ylimmax=105)
            else:
                ylabel = "Root mean squares error"
                std_risk  = np.std(np.asarray(risks),ddof=1,axis=1)
                up_risk = mean_risk +std_risk
                low_risk = mean_risk-std_risk
                risk_plot =  np.asarray((cal_risks,ex_risks,mean_risk,up_risk,low_risk))
                ret = statplot.plotline(np.arange(1,opt.nlvsmax+1),risk_plot,fig_prefix,xlabel,ylabel,['b-','g-','r-','r--','r--'],["Calibration","External Validation","Cross Validation"],title="Latentvariables selection",linewidth=1.5)

    if opt.nlvsoptimize:
        sys.stderr.write("[INFO] start nlvs optimize for %s...\n"%opt.phenotype)
        opt.nlvsmax = min(data.p-1,opt.nlvsmax)
        if check_dimension(data.p,opt.nlvsmax,opt):
            sys.stderr.write("number of variables (p) < nlvsmax")
        risks = np.zeros((opt.nlvsmax,30))
        cal_risks = np.zeros(opt.nlvsmax)
        for i in range(1,opt.nlvsmax+1):
            for j in range(30):
                risk = gwaspls.plscv(data.data,sampleinfo.traits,i,opt.nfold,method=opt.phenotype)
                sys.stderr.write("%s\n"%(str(risk)))
                #print sampleinfo.traits
                risks[i-1,j] = risk
            risk_cal = gwaspls.pls_calibration_risk(data.data,sampleinfo.traits,i,method=opt.phenotype)
            cal_risks[i-1] = risk_cal
            sys.stderr.write("nLVs: %d\n"%i)
        xlabel = "Number of latent variables"
        fig_prefix = opt.outdir+os.path.sep+opt.prefix+".nlvsoptimize"
        xticks_labels = list(map(str,range(1,opt.nlvsmax+1)))
        #ret = statplot.plotyy(np.asarray(svector),risk_plot,numbers,fig_prefix,"Stability",ylabel1,ylabel2,title)
        mean_risk = np.mean(np.asarray(risks),axis=1)
        if opt.phenotype == "category":
            ylabel = "Classification accuracy"
            up_risk = np.max(np.asarray(risks),axis=1)
            low_risk = np.min(np.asarray(risks),axis=1)
            risk_plot =  np.asarray((cal_risks,mean_risk,up_risk,low_risk))
            ret = statplot.plotline(np.arange(1,opt.nlvsmax+1),risk_plot*100,fig_prefix,xlabel,ylabel,['b-','r-','r--','r--'],["calibration","cross validation"],title="Latent variables selection",linewidth=2.0,ylimmax=105)

        else:
            ylabel = "Root mean squares error"
            std_risk  = np.std(np.asarray(risks),ddof=1,axis=1)
            up_risk = mean_risk +std_risk
            low_risk = mean_risk-std_risk
            risk_plot =  np.asarray((cal_risks,mean_risk,up_risk,low_risk))
            ret = statplot.plotline(np.arange(1,opt.nlvsmax+1),risk_plot,fig_prefix,xlabel,ylabel,['b-','r-','r--','r--'],["calibration","cross validation"],title="Latent  variables selection",linewidth=2.0)
    
    if opt.plot:
        result = gwaspls.plsgwas(data.data,sampleinfo.traits,nlvs=opt.ndimentions)
        if opt.ndimentions !=2 and opt.ndimentions != 3:
            sys.stderr.write("Unkown plot dimentions: %dD",opt.ndimentions)
            return 1
        else:
            ret = statplot.plot_Xscore(result.scoreX,sampleinfo.classnums,sampleinfo.uniqclassnum,sampleinfo.uniqcolor,sampleinfo.uniqmarker,sampleinfo.uniqclasslabel,opt.outdir+os.path.sep+opt.prefix+".Xscore","1st latent variable","2nd latent variable","3rd latent variable",dim=opt.ndimentions)
            if opt.phenotype != "category":
                if opt.ndimentions == 3:
                    ret = statplot.plot_XYscore(result.scoreX,sampleinfo.traits,sampleinfo.classnums,sampleinfo.uniqclassnum,sampleinfo.uniqcolor,sampleinfo.uniqmarker,sampleinfo.uniqclasslabel,opt.outdir+os.path.sep+opt.prefix+".XYscore","1st latent variable","2nd latent variable","Traits",3)
                else:
                    ret = statplot.plot_XYscore(result.scoreX,sampleinfo.traits,sampleinfo.classnums,sampleinfo.uniqclassnum,sampleinfo.uniqcolor,sampleinfo.uniqmarker,sampleinfo.uniqclasslabel,opt.outdir+os.path.sep+opt.prefix+".XYscore","1st latent variable","Traits",dim=2)
    if opt.mcuve or opt.viptest:
        calVIPdat(data.data,sampleinfo.traits,opt,data.p)
    #if opt.viptest:
    #    pvalue = do_viptest(data.data,sampleinfo.traits,opt,data.p)
    #    if pvalue is None: return 1
    #    #del data.data
    #    ret_pvalue = viptest_output(fc_rho,pvalue,opt)
    #    if ret_pvalue is None: ret = 1
    if opt.mcuve == opt.mcvip == 1:
        sys.stderr.write("mcuve, mcvip, vip_mcvip should choose one of them")
        return 1
    stability = None
    if opt.mcvip:
        pvalue = do_viptest(data.data,sampleinfo.traits,opt,data.p)
        if pvalue is None: return 1
        ret_pvalue = viptest_output(fc_rho,pvalue,opt)
        del data.data
        del data
        data = MatrixAnno()
        prefix,suffix = utils.parse_filename(os.path.basename(opt.matrix_anno))
        data.parse_matrix_anno(opt.outdir+os.path.sep+prefix+".viptest.anno")
        ret = centring(data.data)
        ret = normalize(data.data)
        if check_dimension(data.p,opt.nlvs,opt):return 1
    if opt.mcuve:
        #print("### do MCUVE")
        stability = gwaspls.plsmcuve(data.data,sampleinfo.traits,opt.mc_nlvs,opt.mcratio,opt.mctimes)
        stability = np.asarray(stability.T)[0]
        stability = np.abs(stability)
        print(stability)
    if opt.mcvip:
        stability = gwaspls.plsmcvip(data.data,sampleinfo.traits,opt.mc_nlvs,opt.mcratio,opt.mctimes)
        stability = np.asarray(stability.T)[0]
    if stability is not None:
        #可以通过粗搜索和精搜索2步, 设计图中图
        best_stat = stability_optimize(data.data,sampleinfo.traits,opt,stability)
        if best_stat is None:
            return 1
        ret = mcmd_output(fc_rho,stability,opt,best_stat)
    return 0


if __name__ == '__main__':
    print(os.path.dirname(os.path.abspath(__file__)))


