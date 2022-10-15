import subprocess as sp
import argparse
import pandas as pd
import numpy as np
import os
import sys

def annotate_gwas_plink2(inbase, outbase, phe, binary=True, gwas_kwds={}, count_kwds={}, freq_kwds={},
                         hwe_kwds={}, dryrun=False, pcut=None):
    '''
    Runs GWAS using --glm function in plink2 and annotates the results file with counts of cases
    and controls, overall maf, and HWE P-value
    
    inbase - str - path to the input pgen
    
    outbase - str - path for output files
    
    phe - str - name of the phenotype
    
    gwas_kwds - dict - contains key value pairs for the --glm plink2 command. Value = None for standalone flag

    count_kwds - dict - contains key value pairs for the --geno-counts plink2 command. Value = None for standalone flag

    freq_kwds - dict - contains key value pairs for the --freq plink2 command. Value = None for standalone flag

    hwe_kwds - dict - contains key value pairs for the --hwe plink2 command. Value = None for standalone flag
    
    dryrun - Boolean - if True, only prints the commands without running them

    Resuts
    --------------------------------------------------------
    Creates an output file <outbase>.<phe>.gwasresults.gwas that should be readable by IGV and returns the final data
    '''
    
    ## Check if variant IDs are unique, if not, error out
    pvar = pd.read_csv(inbase+'.pvar', sep='\t', comment='#', header=None)
    if not pvar[2].is_unique:
        print('Variant IDs must be unique!')
        sys.exit('-1')
    
    
    glmcmd = 'plink2 --glm --pfile {gen} --out {out}'.format(gen=inbase, out=outbase)
    for k in gwas_kwds.keys():
        if gwas_kwds[k] is None:
            glmcmd += ' %s' % k
        else:
            glmcmd += ' %s %s' % (k, gwas_kwds[k])
            
    print(glmcmd+'\n\n')
    
    if not dryrun:
        result = str(sp.check_output(glmcmd, shell=True, stderr=sp.STDOUT))
        print(result)
        
    if binary:
        countctrlcmd = 'plink2 --geno-counts --pfile {gen} --out {out}.{phe}.control --keep-if {phe} == 1'.format(gen=inbase, out=outbase, phe=phe)
        for k in count_kwds.keys():
            if count_kwds[k] is None:
                countctrlcmd += ' %s' % k
            else:
                countctrlcmd += ' %s %s' % (k, count_kwds[k])
            
        print(countctrlcmd+'\n\n')
    
        if not dryrun:
            result = sp.check_output(countctrlcmd, shell=True, stderr=sp.STDOUT)
            print(result)
    
        countcasecmd = 'plink2 --geno-counts --pfile {gen} --out {out}.{phe}.case --keep-if {phe} == 2'.format(gen=inbase, out=outbase, phe=phe)
        for k in count_kwds.keys():
            if count_kwds[k] is None:
                countcasecmd += ' %s' % k
            else:
                countcasecmd += ' %s %s' % (k, count_kwds[k])
            
        print(countcasecmd+'\n\n')
    
        if not dryrun:
            result = sp.check_output(countcasecmd, shell=True, stderr=sp.STDOUT)
            print(result)
    
    freqcmd = 'plink2 --freq --pfile {gen} --out {out}.{phe} --require-pheno {phe}'.format(gen=inbase, out=outbase, phe=phe)
    for k in freq_kwds.keys():
        if freq_kwds[k] is None:
            freqcmd += ' %s' % k
        else:
            freqcmd += ' %s %s' % (k, freq_kwds[k])
            
    print(freqcmd+'\n\n')
    
    if not dryrun:
        result = sp.check_output(freqcmd, shell=True, stderr=sp.STDOUT)
        print(result)
   
    hwecmd = 'plink2 --hardy --pfile {gen} --out {out}.{phe} --require-pheno {phe}'.format(gen=inbase, out=outbase, phe=phe)
    for k in hwe_kwds.keys():
        if hwe_kwds[k] is None:
            hwecmd += ' %s' % k
        else:
            hwecmd += ' %s %s' % (k, hwe_kwds[k])
            
    print(hwecmd+'\n\n')
    
    if not dryrun:
        result = sp.check_output(hwecmd, shell=True, stderr=sp.STDOUT)
        print(result)
        
    if not dryrun:
        if binary == True:
            gwasoutfn = '{out}.{phe}.glm.logistic'.format(out=outbase, phe=phe)
        else:
            gwasoutfn = '{out}.{phe}.glm.linear'.format(out=outbase, phe=phe)
            
        resdata = pd.read_csv(gwasoutfn, sep='\t', dtype={'#CHROM':str})
        resdata = resdata.rename(columns={'#CHROM':'CHR', 'POS':"BP", 'ID':'SNP'}).query('TEST == "ADD"') 
        
        if binary:
            casecount = pd.read_csv('{out}.{phe}.case.gcount'.format(out=outbase, phe=phe), sep='\t', dtype={'#CHROM':str})
            controlcount = pd.read_csv('{out}.{phe}.control.gcount'.format(out=outbase, phe=phe), sep='\t', dtype={'#CHROM':str})
            casecount['casecount'] = casecount.apply(lambda x: '%d | %d | %d' % (x['HOM_REF_CT'], x['HET_REF_ALT_CTS'], x['TWO_ALT_GENO_CTS']), axis=1)
            controlcount['controlcount'] = controlcount.apply(lambda x: '%d | %d | %d' % (x['HOM_REF_CT'], x['HET_REF_ALT_CTS'], x['TWO_ALT_GENO_CTS']), axis=1)
            casecount.rename(columns={'MISSING_CT':'casemissing', 'ID':'SNP'}, inplace=True)
            controlcount.rename(columns={'MISSING_CT':'controlmissing', 'ID':'SNP'}, inplace=True)
            
            X = casecount[['HOM_REF_CT','HET_REF_ALT_CTS','TWO_ALT_GENO_CTS','HAP_REF_CT','HAP_ALT_CTS','casemissing']]
            casecount['casefreq'] = X.div(X.sum(axis=1), axis=0).apply(lambda x: '%.2f | %.2f | %.2f' % (x['HOM_REF_CT'], x['HET_REF_ALT_CTS'], x['TWO_ALT_GENO_CTS']), axis=1)

            X = controlcount[['HOM_REF_CT','HET_REF_ALT_CTS','TWO_ALT_GENO_CTS','HAP_REF_CT','HAP_ALT_CTS','controlmissing']]
            controlcount['controlfreq'] = X.div(X.sum(axis=1), axis=0).apply(lambda x: '%.2f | %.2f | %.2f' % (x['HOM_REF_CT'], x['HET_REF_ALT_CTS'], x['TWO_ALT_GENO_CTS']), axis=1)         
            
            resdata = pd.merge(left=resdata, right=casecount[['SNP', 'casecount', 'casemissing','casefreq']], on='SNP', how='left')
            resdata = pd.merge(left=resdata, right=controlcount[['SNP', 'controlcount', 'controlmissing','controlfreq']], on='SNP', how='left')     
        
        



        freq = pd.read_csv('{out}.{phe}.afreq'.format(out=outbase, phe=phe), sep='\t', dtype={'#CHROM':str})
        ind = freq['ALT_FREQS']>.5
        freq.loc[ind,'ALT_FREQS'] = 1 - freq.loc[ind,'ALT_FREQS']
        freq.rename(columns={'ALT_FREQS':'MAF', 'ID':'SNP'}, inplace=True)
        resdata = pd.merge(left=resdata, right=freq[['SNP', 'MAF']], on='SNP', how='left')     
        
        hwefn = '{out}.{phe}.hardy'.format(out=outbase, phe=phe)
        print(hwefn)
        if os.path.exists(hwefn):
            hwe = pd.read_csv(hwefn, sep='\t', dtype={'#CHROM':str})
            hwe.rename(columns={'P':'HWE_P', 'ID':'SNP'}, inplace=True)
            resdata = pd.merge(left=resdata, right=hwe[['SNP', 'HWE_P']], on='SNP', how='left')     

        resdata.to_csv('{out}.{phe}.gwasresults.gwas'.format(out=outbase, phe=phe), sep='\t', index=False)
        
        if pcut is not None:
            outfn = '{out}.{phe}.gwasresults.Plessthan'.format(out=outbase, phe=phe) + '%.2f.gwas' % pcut
            resdata.loc[resdata['P']<pcut, :].to_csv(outfn, sep='\t', index=False)
            
        return(resdata)
        
        

        