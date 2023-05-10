#!/extdata6/Jinman/opt/python3/bin/python3

## For parsing manifest from TCGA and creating new patient IDs with softlinks to original vcf files

import os, sys, pickle, time, re, copy, subprocess, uuid, math, struct, scipy, itertools
from scipy import stats
from operator import itemgetter
import numpy as np

## region Globals
sTIME_STAMP     = '%s' % (time.ctime().replace(' ', '-').replace(':', '_'))

fPSEUDO_COUNT   = float('1.0e-300')

sBASE_DIR       = '/extdata6/Jinman/07_rbp'
sECLIP_DIR      = '%s/eclipdata'  % sBASE_DIR
sFILEID_DIR     = '%s/00_fileids' % sBASE_DIR
sREF_DIR        = '%s/ref'        % sBASE_DIR

sGENOME         = 'hg19'
sSTUDY          = 'TCGA_LAML'
sSEQDATA        = 'WGS'

list_sCHRIDs    = ['chr%s' % (i+1) for i in range(22)] + ['chrX', 'chrY']
list_sCANCERs   = ['TCGA_LAML', #'TCGA_LIHC',
                   'TCGA_BRCA', 'TCGA_LUAD',
                   'TCGA_UCEC', 'TCGA_OV',
                   'TCGA_COAD', 'TCGA_CESC',
                   'TCGA_SARC', 'TCGA_ESCA']

dict_sCELLS     = {'TCGA_LIHC': 'HepG2','TCGA_LAML': 'K562',
                   'TCGA_BRCA': 'Both', 'TCGA_LUAD': 'Both',
                   'TCGA_UCEC': 'Both', 'TCGA_OV'  : 'Both',
                   'TCGA_COAD': 'Both', 'TCGA_CESC': 'Both',
                   'TCGA_SARC': 'Both', 'TCGA_ESCA': 'Both'
                   }

dict_sSAMPLECNT = {'TCGA_LIHC': 378, 'TCGA_LAML': 182,
                   'TCGA_BRCA': 510, 'TCGA_LUAD': 479,
                   'TCGA_UCEC': 528, 'TCGA_OV'  : 441,
                   'TCGA_COAD': 188, 'TCGA_CESC': 228,
                   'TCGA_SARC': 173, 'TCGA_ESCA': 305
                   }

dict_sCAPTURE   = {'TCGA_LIHC': '%s/VCRome21.bed'                % sREF_DIR,
                   'TCGA_LAML': '%s/Agilent_SureSelect_50Mb.bed' % sREF_DIR,
                   'TCGA_LUAD': '%s/Agilent_SureSelect_V4.bed'   % sREF_DIR,
                   'TCGA_BRCA': '%s/SeqCapExomeV3.bed'           % sREF_DIR,
                   'TCGA_UCEC': '', 'TCGA_OV'  : '',
                   'TCGA_COAD': '', 'TCGA_CESC': '',
                   'TCGA_SARC': '', 'TCGA_ESCA': ''
                   }


sCELLS          = dict_sCELLS[sSTUDY]
nSAMPLE_CNT     = dict_sSAMPLECNT[sSTUDY]

sCAPTUREFILE    = dict_sCAPTURE[sSTUDY]
sREFINDEX       = '/extdata6/Jinman/reference_genome/hg19/hg19.fa.fai'
sREFSEQFILE     = '/extdata6/Jinman/reference_genome/091613_4WT.txt'
sANNOVAR        = '/extdata6/Jinman/util/annovar/table_annovar.pl'

s3UTR_FILE      = { 'hg19'      : '/extdata6/Jinman/07_rbp/ref/3UTRseq_hg19.txt'
                    }

dict_sGENOME    = { 'hg19'      :'/extdata6/Jinman/reference_genome/hg19/hg19.fa',
                    }

dict_sMIRBASE   = {'hg19': '/extdata6/Jinman/reference_genome/mirbase/hsa.liftover.gff3'
                    }

dict_sMIRFILE   = {'HepG2':'%s/top20_miRNAs_HepG2.txt' % sREF_DIR,
                   'K562':'%s/top20_miRNAs_K562.txt'   % sREF_DIR}

dict_sGENIC     = {'exonic': 'exonic', 'UTR3': '3UTR', 'UTR5': '5UTR', 'intronic': 'intronic',
                   'splicing': 'splicing', 'ncRNA_splicing': 'splicing', 'downstream': 'downstream',
                   'ncRNA_exonic': 'intergenic', 'intergenic': 'intergenic', 'ncRNA_UTR5': '5UTR',
                   'ncRNA_UTR3': '3UTR', 'upstream': 'upstream', 'ncRNA_intronic': 'intergenic'}

dict_sGENIC_alt = {'exonic': 'exonic', 'UTR3': '3UTR', 'UTR5': '5UTR', 'intronic': 'intronic',
                   'splicing': 'other', 'ncRNA_splicing': 'other', 'downstream': 'intergenic',
                   'ncRNA_exonic': 'other', 'intergenic': 'intergenic', 'ncRNA_UTR5': '5UTR',
                   'ncRNA_UTR3': '3UTR', 'upstream': 'other', 'ncRNA_intronic': 'other', '5UTR': '5UTR',
                   '3UTR':'3UTR', 'other':'other'}

list_sAllST     = ['8mer', '7mer-m8', '7mer-A1', '6mer',
                   '6mer-A1', 'o7mer', 'o6mer',
                   'CDNST_1', 'CDNST_2','CDNST_3', 'CDNST_4']

dict_WO_ST      = {'8mer':    ['7mer-m8', '7mer-A1', '6mer'],
                 '7mer-m8': ['8mer', '7mer-A1', '6mer'],
                 '7mer-A1': ['8mer', '7mer-m8', '6mer'],
                 '6mer':    ['8mer', '7mer-m8', '7mer-A1'],

                 '6mer-A1': ['8mer', '7mer-m8', '7mer-A1', '6mer', 'o7mer', 'o6mer'],
                 'o7mer':   ['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o6mer'],
                 'o6mer':   ['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o7mer'],

                 'pivot':   ['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o7mer',
                             'o6mer', 'centered1', 'centered2', 'CDNST_1', 'CDNST_2', 'CDNST_3', 'CDNST_4', 'CDNST_5'],
                 'centered1':['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o7mer',
                             'o6mer', 'pivot', 'centered2','CDNST_1', 'CDNST_2', 'CDNST_3', 'CDNST_4', 'CDNST_5'],
                 'centered2':['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o7mer',
                             'o6mer', 'pivot', 'centered1', 'CDNST_1', 'CDNST_2', 'CDNST_3', 'CDNST_4', 'CDNST_5'],

                 'CDNST_1': ['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o7mer', 'o6mer',
                             'CDNST_2', 'CDNST_3', 'CDNST_4'],
                 'CDNST_2': ['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o7mer', 'o6mer',
                             'CDNST_1', 'CDNST_3', 'CDNST_4'],
                 'CDNST_3': ['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o7mer', 'o6mer',
                             'CDNST_1', 'CDNST_2', 'CDNST_4'],
                 'CDNST_4': ['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o7mer', 'o6mer',
                             'CDNST_1', 'CDNST_2', 'CDNST_3'],
                 #'CDNST_5': ['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o7mer', 'o6mer',
                  #           'CDNST_1', 'CDNST_2', 'CDNST_3', 'CDNST_4'],
                 'NC':      ['8mer', '7mer-m8', '7mer-A1', '6mer', '6mer-A1', 'o7mer', 'o6mer',
                             'CDNST_1', 'CDNST_2', 'CDNST_3', 'CDNST_4']
               }


list_sCST       = ['8mer', '7mer-m8', '7mer-A1', '6mer']
list_s78MER     = ['8mer', '7mer-m8', '7mer-A1']
list_sNST       = ['6mer-A1', 'o7mer', 'o6mer']
list_sCDNST     = ['CDNST_1', 'CDNST_2', 'CDNST_3', 'CDNST_4']

'''V-S Check Thresholds'''
nREFSEQ_GENECNT             = 18625

n3UTR_FILE_COLUMN_SIZE      = 3
nMAX_3UTR_SIZE              = 1000    ## Max 3UTR Size for ALL
nGENOME_SIZE                = 3000000000
nGFF3_COLUMN_SIZE           = 9

'''Minimum number of variants in cluster '''
nMIN_VARIANTS               = 2

nMAX_3UTR_SIZE              = 1000    ## Max 3UTR Size for ALL
nMIN_NO_SITE_SIZE           = 200     ## MIN 3UTR Size for just No_Site UTRS
fTOP_EXPR                   = 1       ## Top Expressed


## endregion

## region Util Functions

def get_clipdata_list (sStudy=sSTUDY, bFull=0):

    '''Load CLIP data'''
    sClipMetaFile  = '%s/metadata.tsv'  % sECLIP_DIR
    list_cClipData = eCLIP_parse_eclip_meta (sClipMetaFile)

    if bFull: return list_cClipData
    else:
        if sStudy in ['TCGA_LIHC', 'TCGA_LAML']:
            list_cClipData_filt = [cClipData for cClipData in list_cClipData if cClipData.sCellLine == sCELLS] #HepG2 (contains K562)

        else: # Intersection of RBPs

            dict_sRBPID = {}

            for cClipData in list_cClipData:
                if cClipData.sGeneSym not in dict_sRBPID:
                    dict_sRBPID[cClipData.sGeneSym] = []
                dict_sRBPID[cClipData.sGeneSym].append(cClipData)
            #loop END: cClipData

            list_cClipData_filt = []
            for sRBPID in dict_sRBPID:
                if len(dict_sRBPID[sRBPID]) == 4:
                    list_cClipData_filt += dict_sRBPID[sRBPID]
            #loop END: sRBPID
        #if END:
        return list_cClipData_filt
    #if END:
#def END: get_clipdata_list

def get_clipdata_list_v2 (sFull=0):

    ## Returns Dictionary of CancerTypes and list of RBPs

    '''Load CLIP data'''
    sClipMetaFile  = '%s/metadata.tsv'  % sECLIP_DIR
    list_cClipData = eCLIP_parse_eclip_meta (sClipMetaFile)

    if sFull: # Full List of RBPs
        dict_sOutput = {}
        dict_sRBPID  = {}

        for cClipData in list_cClipData:
            if cClipData.sGeneSym not in dict_sRBPID:
                dict_sRBPID[cClipData.sGeneSym] = []
            dict_sRBPID[cClipData.sGeneSym].append(cClipData)
        #loop END: cClipData

        for sCancer in list_sCANCERs:
            if sCancer not in dict_sOutput:
                dict_sOutput[sCancer] = []
            dict_sOutput[sCancer] = list(dict_sRBPID.keys())

    else: # RBPs according to tissue or intersection
        dict_sOutput        = {}

        for sCancer in list_sCANCERs:
            list_cClipData_filt = []

            if sCancer not in dict_sOutput:
                dict_sOutput[sCancer] = {}

            if sCancer == 'TCGA_LIHC':
                list_cClipData_filt = [cClipData for cClipData in list_cClipData if cClipData.sCellLine == 'HepG2']

            elif sCancer == 'TCGA_LAML':
                list_cClipData_filt = [cClipData for cClipData in list_cClipData if cClipData.sCellLine == 'K562']

            else:
                dict_sRBPID = {}
                for cClipData in list_cClipData:
                    if cClipData.sGeneSym not in dict_sRBPID:
                        dict_sRBPID[cClipData.sGeneSym] = []
                    dict_sRBPID[cClipData.sGeneSym].append(cClipData)
                #loop END: cClipData

                for sRBPID in dict_sRBPID:
                    if len(dict_sRBPID[sRBPID]) == 4:
                        list_cClipData_filt += dict_sRBPID[sRBPID]
                #loop END: sRBPID
            #if END:

            for cClipData in list_cClipData_filt:
                if cClipData.sGeneSym not in dict_sOutput[sCancer]:
                    dict_sOutput[sCancer][cClipData.sGeneSym] = []
                dict_sOutput[sCancer][cClipData.sGeneSym].append(cClipData.sAccID)
           #loop END: cClipData
        # loop END: sCancer
    #if END:
    return dict_sOutput
#def END: get_clipdata_list

def get_chrom_sizes():

    dict_sChrSize = {}

    for sReadLine in open(sREFINDEX, 'r'):
        list_sCol   = sReadLine.split('\t')
        sChrID      = list_sCol[0].replace('chr', '')
        nChrSize    = int(list_sCol[1])

        if sChrID not in dict_sChrSize:
            dict_sChrSize[sChrID] = 0
        dict_sChrSize[sChrID] = nChrSize
    # loop END: sReadLine

    return dict_sChrSize
#def END: get_chrom_sizes

def set_chr_range(sChrID, nWinSize):

    nChrSize     =  get_chrom_sizes()[sChrID[3:]]
    list_nWindow = []
    nRange       = int(nChrSize/nWinSize) + 1

    for i in range(nRange):

        x = i * nWinSize + 1
        y = (i+1) * nWinSize

        if y > nChrSize: y = nChrSize

        list_nWindow.append('%s-%s' % (int(x), int(y)))
    #loop END: i
    return list_nWindow
#def END: set_chr_range

def get_patientIDs (sStudy, sSeq):

    sInFile         = '%s/%s_%s.txt' % (sFILEID_DIR, sStudy, sSeq)
    InFile          = open(sInFile, 'r')
    list_sFileIDs   = []

    for sReadLine in InFile:
        list_sFileIDs.append (sReadLine.strip('\n'))
    #loop END:

    return list_sFileIDs
#def END: get_file_ids

def histogram (list_sInput):

    nBins        = int(max(list_sInput)) + 1
    ## MAF Histogram ##
    dict_fMAF = {i:[] for i in range(nBins)}

    fBinRange = 5

    for nGroup, lCnt in itertools.groupby(list_sInput, key=lambda n: n//fBinRange):
        dict_fMAF[nGroup] = dict_fMAF[nGroup] + list(lCnt)
    #loop END: nGroup, lCnt

    for nGroup in dict_fMAF:
        print(nGroup, len(dict_fMAF[nGroup]))
    #loop END: nGroup
#def END: output_for_histogram

def reverse_complement(sSeq):
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '.':'.'}
    list_sSeq       = list(sSeq.upper()) # Turns the sequence in to a gigantic list
    list_sSeq       = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1] # Empty string, join the list of bases, [start : end : backwards]
#def END: reverse_complement

def process_chromo_ID(sChrID):
    sChrID = sChrID.upper()

    if '_' not in sChrID:  # Only with those chromosomes that are distinctly 1-22 or X and Y
        sChrID = sChrID.strip('CHR')

        if sChrID == 'X':
            sChrID = '23'
            return int(sChrID)

        elif sChrID == 'Y':
            sChrID = '24'
            return int(sChrID)
        else:
            return int(sChrID)
    else:
        sChrID = "0"  # Flagged as misc chromosome
        return int(sChrID)
# def END: process_chromo_ID

def copy_temp_core_script (sWorkDir):
    os.makedirs('%s/temp' % sWorkDir, exist_ok=True)
    os.system('cp %s/A_main.py %s/temp/A_main_%s.py'
              % (sWorkDir, sWorkDir, sTIME_STAMP))
    return '%s/temp/A_main_%s.py' % (sWorkDir, sTIME_STAMP)
#def END: copy_temp_core_script

def print_start_msg (sTaskName, sFileDir=''):
    print('START...%s %s %s' % (sTaskName, sFileDir, time.ctime()))
#def END: print_start_msg

def print_done_msg (sTaskName, sFileDir=''):
    print('DONE...%s %s %s' % (sTaskName, sFileDir, time.ctime()))
#def END: print_done_msg

def process_chromo_ID(sChrID):
    sChrID = sChrID.upper().replace('CHR', '')
    if sChrID == 'X':   return 23
    elif sChrID == 'Y': return 24
    elif sChrID == 'M': return 25
    else:               return int(sChrID)
#def END: process_chromo_ID

def load_reference_gencode (sInFile):

    InFile       = open(sInFile, 'r')
    dict_sOutput = {}
    for i,sReadLine in enumerate(InFile):
        if sReadLine.startswith('##'): continue #Skip headers

        ''' File Format
        0 sChrID          chr1
        1 sAnnotation     HAVANA
        2 sFeature        gene
        3 nStartPos       11869   *1-based
        4 nEndPos         14409
        5 sScore          .       *Not used
        6 sStrand         +
        7 nGenomicPhase   .
        8 sAddInfo        gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1";  level 2; havana_gene "OTTHUMG00000000961.2";
        '''
        list_sColumn    = sReadLine.strip('\n').split('\t')

        sChrID          = list_sColumn[0]
        sAnnotation     = list_sColumn[1]
        sFeature        = list_sColumn[2]

        if sFeature != 'gene': continue # Skip other redundant info (exon, CDS, transcript, UTR...)

        nStartPos       = int(list_sColumn[3])
        nEndPos         = int(list_sColumn[4])
        sScore          = list_sColumn[5]
        sStrand         = list_sColumn[6]
        nGenomicPhase   = list_sColumn[7]
        sGeneInfo       = list_sColumn[8].replace('"', '')

        dict_sGeneInfo  = dict([list(filter(None, sInfo.split(' '))) for sInfo in sGeneInfo.split(';') if list(filter(None, sInfo.split(' ')))])

        sGenesym        = dict_sGeneInfo['gene_name']
        sGeneID         = dict_sGeneInfo['gene_id']

        if sGeneID not in dict_sOutput:
            dict_sOutput[sGeneID] = ''
        dict_sOutput[sGeneID] = sGenesym

    #loop END: sReadLine

    #VS-Check
    if not dict_sOutput:
        sys.exit('Invalid Dict : load_reference_gencode : dict_sOutput size= %d' % len(dict_sOutput))

    return dict_sOutput
#def END: load_reference_gencode

def count_site_types (sTempSeq, sSiteSeq):

    list_sSeq = list(sTempSeq)
    list_nPos = []
    nMotifCnt = 0

    for sReIndex in re.finditer(sSiteSeq, sTempSeq):
        nIndexStart = sReIndex.start()
        nIndexEnd   = sReIndex.end()

        list_nPos.append([nIndexStart, nIndexEnd-1])

        nMotifCnt   += 1

        for i in range(nIndexStart, nIndexEnd):
            list_sSeq[i] = '*'
        #loop END: i
    #loop END: sReIndex

    sNewSeq = ''.join(list_sSeq)

    #V-S Check: String Length
    if len(sTempSeq) != len(sNewSeq):
        sys.exit('ERROR: sTempSeq Size= %d\tsNewSeq Size= %d' % (len(sTempSeq), len(sNewSeq)))

    return sNewSeq, nMotifCnt, list_nPos
#def END: count_site_types

def determine_site_types (sMirSeq, sNegControl = '1'):

    #V-S Check: Valid Sequence
    if 'N' in sMirSeq or len(sMirSeq) == 0:
        sys.exit('ERROR: Invalid MirSeq %s' % sMirSeq)

    dict_sSTs            = {}
    dict_sWobble         = {'A':'A', 'C':'C', 'G':'A', 'T':'C'}


    ## CST ##
    dict_sSTs['8mer']    = reverse_complement(sMirSeq[2-1:8].replace('U','T')).replace('T','U')+'A'  # 2-8 match + A
    dict_sSTs['7mer-m8'] = reverse_complement(sMirSeq[2-1:8].replace('U','T')).replace('T','U')      # 2-8 match
    dict_sSTs['7mer-A1'] = reverse_complement(sMirSeq[2-1:7].replace('U','T')).replace('T','U')+'A'  # 2-7 match + A
    dict_sSTs['6mer']    = reverse_complement(sMirSeq[2-1:7].replace('U','T')).replace('T','U')      # 2-7 match

    ## NST ##
    dict_sSTs['6mer-A1'] = reverse_complement(sMirSeq[2-1:6].replace('U','T')).replace('T','U')+'A'  # 2-6 match + A
    dict_sSTs['o7mer']   = reverse_complement(sMirSeq[3-1:9].replace('U','T')).replace('T','U')      # 3-9 match
    dict_sSTs['o6mer']   = reverse_complement(sMirSeq[3-1:8].replace('U','T')).replace('T','U')      # 3-8 match
    dict_sSTs['o5mer']   = reverse_complement(sMirSeq[3-1:7].replace('U','T')).replace('T','U')      # 3-7 match

    ## PRNST ##
    dict_sSTs['pivot']     = reverse_complement(sMirSeq[2-1:4].replace('U','T') + '.' + sMirSeq[6-1:7].replace('U','T')).replace('T','U')
    dict_sSTs['centered1'] = reverse_complement(sMirSeq[4-1:14].replace('U','T')).replace('T','U')
    dict_sSTs['centered2'] = reverse_complement(sMirSeq[5-1:15].replace('U','T')).replace('T','U')

    ## CDNST ##
    dict_sSTs['5mer']    = reverse_complement(sMirSeq[2-1:6].replace('U','T')).replace('T','U')      # 2-6 match

    dict_sSTs['CDNST_1'] = reverse_complement(sMirSeq[2-1:6].replace('U','T')).replace('T','U')      # 2-6 match
    dict_sSTs['CDNST_2'] = reverse_complement(sMirSeq[6-1:7].replace('U','T')).replace('T','U') + '.' +\
                           reverse_complement(sMirSeq[2-1:4].replace('U','T')).replace('T','U') +'A'

    dict_sSTs['CDNST_3'] = reverse_complement(sMirSeq[7-1:9].replace('U','T')).replace('T','U') + '.' +\
                           reverse_complement(sMirSeq[4].replace('U','T')).replace('T','U') + '.' +\
                           reverse_complement(sMirSeq[2].replace('U','T')).replace('T','U')

    dict_sSTs['CDNST_4'] = reverse_complement(sMirSeq[7].replace('U','T')).replace('T','U') + '...' +\
                           reverse_complement(sMirSeq[2-1:4].replace('U','T')).replace('T','U') + 'A'

    dict_sSTs['CDNST_5'] = reverse_complement(sMirSeq[8].replace('U','T')).replace('T','U') + '..' +\
                           reverse_complement(sMirSeq[2-1:5].replace('U','T')).replace('T','U')

    ## Negative Control ##
    sSite_NegCon1        = reverse_complement(sMirSeq[10-1:16].replace('U','T'))    # 10-15 match
    sSite_NegCon2        = reverse_complement(sMirSeq[6-1:12].replace('U','T'))     # 6-11 match
    sSite_NegCon3        = reverse_complement(sMirSeq[5-1:11].replace('U','T'))     # 5-10 match

    # Sample Mir:         _T_A_C_A_T_A_C_T_T_C_T_T_T_A_C_A_T_T_C_C_A
    # Negative Control 1:            A_C_T_._C_._T
    # Negative Control 3:              ._T_T_._._C_A_
    # Negative Control 3:                      T_T_._._._A_T

    # Negative Control 1:  5'X_X_X_X_X_X_X_X_O_O_O_W_O_W_O_X_3'
    sSite_NegCon1        = sSite_NegCon1[:3] + dict_sWobble[sSite_NegCon1[3:4]] +  sSite_NegCon1[4:5] + dict_sWobble[sSite_NegCon1[5:6]] + sSite_NegCon1[6:]

    # Negative Control 2:  5'X_X_X_X_X_X_O_W_O_W_O_O_O_X_X_X_3'
    sSite_NegCon2        = sSite_NegCon2[:1] + dict_sWobble[sSite_NegCon2[1:2]] + sSite_NegCon2[2:3] + dict_sWobble[sSite_NegCon2[3:4]] + sSite_NegCon2[4:]

    # Negative Control 3:  5'X_X_X_O_O_W_W_W_O_O_X_X_X_X_X_X_3'
    sSite_NegCon3        = sSite_NegCon3[:2] + '%s%s%s' % (dict_sWobble[sSite_NegCon3[2]], dict_sWobble[sSite_NegCon3[3]], dict_sWobble[sSite_NegCon3[4]] ) + sSite_NegCon3[5:]

    if sNegControl   == '1': dict_sSTs['NC']     = sSite_NegCon1.replace('T','U')
    elif sNegControl == '2': dict_sSTs['NC']     = sSite_NegCon2.replace('T','U')
    elif sNegControl == '3': dict_sSTs['NC']     = sSite_NegCon3.replace('T','U')
    else: sys.exit('ERROR: Invalid NegControl Designation %s' % sNegControl)


    ## w5mer # 5'._._._._W_O_O_O_O_._3'

    sMirSeq_w5mer       = sMirSeq[5-1:9] # 5-9 match, wobble on 1st position

    # w5mer
    if sMirSeq_w5mer[0] == 'G':
        sMirSeq_w5mer = 'A' + sMirSeq_w5mer[1:]
        dict_sSTs['w5mer'] = reverse_complement(sMirSeq_w5mer.replace('U','T')).replace('T','U')


    elif sMirSeq_w5mer[0] == 'T':
        sMirSeq_w5mer = 'C' + sMirSeq_w5mer[1:]
        dict_sSTs['w5mer'] = reverse_complement(sMirSeq_w5mer.replace('U','T')).replace('T','U')

    else:
        dict_sSTs['w5mer'] = 'NULL'
    #if END: sMirSeq_w5mer

    ## w9mer   - Survey for potential noise         # 5'._._._._._W_._O_O_O_._._W_O_3'
    sMirSeq_w9mer       = sMirSeq[6-1:14] # 5-9 match, wobble on 1st position

    if sMirSeq_w9mer[0] in ['G','T'] and sMirSeq_w9mer[7] in ['G','T']:

        sMirSeq_w9mer = dict_sWobble[sMirSeq_w9mer[0]] + '.' + sMirSeq_w9mer[2:5] + '..' + dict_sWobble[sMirSeq_w9mer[7]] + sMirSeq_w9mer[8]

        dict_sSTs['w9mer'] = reverse_complement(sMirSeq_w9mer.replace('U','T'))

    else:
        dict_sSTs['w9mer'] = 'NULL'

    ## CD_NST(O_._O_O_O_A)  * Old
    dict_sSTs['SM-6mer-A1']  = reverse_complement(sMirSeq[5].replace('U','T')).replace('T','U') + '.' \
                              + reverse_complement(sMirSeq[2-1:4].replace('U','T')).replace('T','U') +'A'

    return dict_sSTs
#def END: determine_site_types

def parse_capture_bedfile (sInFile):
    print('Loading capture file', sInFile)
    InFile        = open(sInFile, 'r')
    dict_sGeneSym = {}
    dict_sNMID    = {}

    for i, sReadLine in enumerate(InFile):
        ''' File Format
        0	sChrID                  chr1
        1	nStartPos               30317 (0-based)
        2	nEndPos                 30527 (1-based)
        3	dict_sGeneInfo          ensembl_gene_id=ENSG00000243485;vega_gene_id=RP11-34P133;...
        '''
        list_sColumn    = sReadLine.strip('\n').split('\t')
        sChrID          = list_sColumn[0]
        nStartPos       = int(list_sColumn[1])
        nEndPos         = int(list_sColumn[2])

        # Different capture file formats
        if sSTUDY in  ['TCGA_LIHC', 'TCGA_BRCA']:
            list_sGeneInfo  = [sInfo.split('=') for sInfo in list_sColumn[3].split(';') if sInfo] if list_sColumn[3] != '.' else []

        elif sSTUDY in ['TCGA_LAML', 'TCGA_LUAD']:
            list_sGeneInfo = [sInfo.split('|') for sInfo in list_sColumn[3].split(',') if sInfo] if list_sColumn[3] not in ['.', '-'] else []

        list_sGeneInfo  = [sInfo for sInfo in list_sGeneInfo if len(sInfo) == 2]

        dict_sGeneInfo  = {}
        for sLabel, sID in list_sGeneInfo:
            if sLabel not in dict_sGeneInfo:
                dict_sGeneInfo[sLabel] = []
            dict_sGeneInfo[sLabel].append(sID.upper())

        if sSTUDY in  ['TCGA_LIHC', 'TCGA_BRCA']:
            try:
                list_sGeneSym = dict_sGeneInfo['refGene_gn']
                list_sNMID    = dict_sGeneInfo['refGene_tx']

            except KeyError: continue
        elif sSTUDY in ['TCGA_LAML', 'TCGA_LUAD']:
            try: dict_sGeneInfo['ref']
            except KeyError: continue

            list_sGeneSym     = [sID for sID in dict_sGeneInfo['ref'] if not sID.startswith('NM_') and not sID.startswith('NR_')]
            list_sNMID        = [sID for sID in dict_sGeneInfo['ref'] if sID.startswith('NM_')]

        sPosKey = '%s,%s' % (nStartPos, nEndPos)

        for sGeneSym in list_sGeneSym:
            if sGeneSym not in dict_sGeneSym:
                dict_sGeneSym[sGeneSym] = {}

            if sChrID not in dict_sGeneSym[sGeneSym]:
                dict_sGeneSym[sGeneSym][sChrID] = []
            dict_sGeneSym[sGeneSym][sChrID].append(sPosKey)

        for sNMID in list_sNMID:
            if sNMID not in dict_sNMID:
                dict_sNMID[sNMID] = {}

            if sChrID not in dict_sNMID[sNMID]:
                dict_sNMID[sNMID][sChrID] = []
            dict_sNMID[sNMID][sChrID].append(sPosKey)
        #loop END: sGeneSym, sNMID
    # loop END: i, sReadLine
    InFile.close()

    # VS Check
    if not dict_sGeneSym:
        sys.exit('Invalid Dictionary : parse_capture_bedfile : dict_sGeneSym size= %d' % len(dict_sGeneSym))

    # VS Check
    if not dict_sNMID:
        sys.exit('Invalid Dictionary : parse_capture_bedfile : dict_sNMID size= %d' % len(dict_sNMID))


    return dict_sGeneSym, dict_sNMID
#def END: parse_capture_bedfile


## region Q-value Correction
def get_qval_list(list_fPvalues):

    list_fQvalues = []
    list_fPvalues = np.array(list_fPvalues)
    nListLen      = float(list_fPvalues.shape[0])
    list_nValues  = [[float(fPval), nID] for fPval, nID in list_fPvalues]
    list_nValues.sort()
    list_nValues.reverse()

    for i, (fPval, nID) in enumerate(list_nValues):
        nRank = nListLen - i
        fQval =(nListLen / nRank) * fPval

        if fQval >= 1: fQval = 1  # No value exceeding 1

        list_fQvalues.append([int(nID), fQval])
    #loop END: i, [fPval, nIndex, sID]

    return list_fQvalues
#def END: calculate_qvalues

## endregion

## region Add Colors

def red(string, e=0):     return '\033[%s31m%s\033[m'%('' if e == 0 else '1;', string)
def green(string, e=0):   return '\033[%s32m%s\033[m'%('' if e == 0 else '1;', string)
def yellow(string, e=0):  return '\033[%s33m%s\033[m'%('' if e == 0 else '1;', string)
def blue(string, e=0):    return '\033[%s34m%s\033[m'%('' if e == 0 else '1;', string)
def magenta(string, e=0): return '\033[%s35m%s\033[m'%('' if e == 0 else '1;', string)
def cyan(string, e=0):    return '\033[%s36m%s\033[m'%('' if e == 0 else '1;', string)
def white(string, e=0):   return '\033[%s37m%s\033[m'%('' if e == 0 else '1;', string)

def get_color (cMir, sResidue):
    if sResidue == cMir.sAltNuc: sResidue = red(sResidue,1)
    elif sResidue == cMir.sRefNuc: sResidue = green(sResidue,1)
    else: sResidue = blue(sResidue,1)
    return sResidue
#def END: get_color
## endregion

## endregion

## region Classes

## region class cFasta
re_nonchr = re.compile('[^a-zA-Z]')
class cFasta:
    def __init__(self, sRefFile):

        #V-S Check: File Existence
        if not os.path.isfile(sRefFile):
           sys.exit('(): File does not exist')

        self.InFile     = open(sRefFile, 'r')
        self.sChrIDList = []
        self.nChromLen  = []
        self.nSeekPos   = []
        self.nLen1      = []
        self.nLen2      = []

        #V-S Check: File Existence
        if not os.path.isfile('%s.fai'%sRefFile):
           sys.exit('.fai file does not exist')

        InFile = open('%s.fai' % sRefFile, 'r')
        for sLine in InFile:
           list_sColumn = sLine.strip('\n').split() # Goes backwards, -1 skips the new line character

           self.sChrIDList.append  (list_sColumn[0])
           self.nChromLen.append   (int(list_sColumn[1]))
           self.nSeekPos.append    (int(list_sColumn[2]))
           self.nLen1.append       (int(list_sColumn[3]))
           self.nLen2.append       (int(list_sColumn[4]))
        #loop END: sLINE
        InFile.close()
        self.sType = []
    #def END: __init_

    def fetch(self, sChrom, sGeneSym = None, nFrom = None, nTo = None, sStrand = '+'):
        assert sChrom in self.sChrIDList, sChrom
        nChrom = self.sChrIDList.index(sChrom)

        if sGeneSym == None: sGeneSym = ''
        if nFrom == None: nFrom = 0
        if nTo   == None: nTo = self.nChromLen[nChrom]

        try: assert(0 <= nFrom) and(nFrom < nTo) and(nTo <= self.nChromLen[nChrom])
        except AssertionError:
            print('cFasta fetch assertion error',sGeneSym, nFrom, nTo)
            sys.exit()

        nBlank = self.nLen2[nChrom] - self.nLen1[nChrom]

        nFrom  = int(nFrom +(nFrom / self.nLen1[nChrom]) * nBlank) # Start Fetch Position

        nTo    = int(nTo   +(nTo   / self.nLen1[nChrom]) * nBlank) # End Fetch Position

        self.InFile.seek(self.nSeekPos[nChrom] + nFrom)            # Get Sequence

        sFetchedSeq = re.sub(re_nonchr, '', self.InFile.read(nTo - nFrom)).upper()

        if   sStrand == '+':
            return sFetchedSeq

        elif sStrand == '-':
            return reverse_complement(sFetchedSeq)

        else:
            sys.exit('Error: moudle seq: invalid strand')
        #if END: sStrand
    #def END: fetch
#class END: Fasta
## endregion

## region class eCLIP
class eCLIP:
    def __init__(self):
        eCLIP.sAccID = ''
        eCLIP.sFileFormat = ''
        eCLIP.sOutputType = ''
        eCLIP.sAssay = ''

        # CellLine Info
        eCLIP.sCellLineID = ''
        eCLIP.sCellLine = ''
        eCLIP.sCellType = ''
        eCLIP.sCellStage = ''
        eCLIP.sGender = ''
        eCLIP.sOrganism = ''
        eCLIP.sTreatment = ''
        eCLIP.sSubFraction = ''
        eCLIP.sPhase = ''

        # Target Protein
        eCLIP.sGeneSym = ''
        eCLIP.sAntibody = ''
        eCLIP.sLibrary = ''
        eCLIP.sLibDepleted = ''
        eCLIP.sExtractMeth = ''
        eCLIP.sLysisMeth = ''
        eCLIP.sXlinkMeth = ''
        eCLIP.sReleaseDate = ''
        eCLIP.sProject = ''
        eCLIP.sProteinConc = ''
        eCLIP.sLibfrag = ''
        eCLIP.sLibSize = ''
        eCLIP.sAge = ''
        eCLIP.sRepCnt = ''
        eCLIP.sTechCnt = ''
        eCLIP.sReadLen = ''
        eCLIP.sRunType = ''
        eCLIP.sPairedEnd = ''
        eCLIP.sPairedWith = ''
        eCLIP.sDerivedFrom = ''
        eCLIP.sSize = ''
        eCLIP.sLabName = ''
        eCLIP.sCheckSum = ''
        eCLIP.sDownURL = ''
        eCLIP.sAssembly = ''
        eCLIP.sPlatform = ''
        # def END: __init__

def eCLIP_parse_eclip_meta(sInFile):

    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for i, sReadLine in enumerate(InFile):
        ''' File Format
        0	File accession	                    ENCFF990XRM
        1	File format	                        bed narrowPeak
        2	Output type	                        peaks
        3	Experiment accession	            ENCSR441YTO
        4	Assay	                            iCLIP
        5	Biosample term id	                EFO:0002067
        6	Biosample term name	                K562
        7	Biosample type	                    immortalized cell line
        8	Biosample life stage	            adult
        9	Biosample sex	                    female
        10	Biosample organism	                Homo sapiens
        11	Biosample treatments
        12	Biosample subcellular fraction term
        13	Biosample phase
        14	Biosample synchronization stage
        15	Experiment target	                TIAL1-human
        16	Antibody accession	                ENCAB050ZSG
        17	Library made from	                RNA
        18	Library depleted in
        19	Library extraction method	        see document
        20	Library lysis method	            see document
        21	Library crosslinking method	        ultraviolet irradiation
        22	Experiment date released	        9/3/2014
        23	Project	                            ENCODE
        24	RBNS protein concentration
        25	Library fragmentation method	    see document
        26	Library size range	                177-207
        27	Biosample Age	                    53 year
        28	Biological replicate(s)	            2
        29	Technical replicate	                1
        30	Read length
        31  Mapped read legnth
        32	Run type
        33	Paired end
        34	Paired with
        35	Derived from	                    ENCFF002DNO
        36	Size	                            10868
        37	Lab	                                ENCODE Consortium
        38	md5sum	                            91efcfbb960ed08917ab50141f39886b
        39	File download URL	                https://www.encodeproject.org/files/ENCFF990XRM/@@download/ENCFF990XRM.bed.gz
        40	Assembly	                        hg19
        41	Platform
        '''

        #if sReadLine.startswith('File'):
            #list_sHeaders = sReadLine.strip('\n').split('\t')
            #continue

        list_sColumn = sReadLine.strip('\n').split('\t')

        #for sHeader,sExample in zip(list_sHeaders, list_sColumn):
            #print('%s,%s' % (sHeader, sExample))
        #if i == 1:break
        cClipData               = eCLIP()

        cClipData.sAccID        = list_sColumn[0]
        cClipData.sFormat       = list_sColumn[1]
        cClipData.sAssay        = list_sColumn[4]

        # CellLine Info
        cClipData.sCellLineID   = list_sColumn[5]
        cClipData.sCellLine     = list_sColumn[6]

        # Target Protein
        cClipData.sGeneSym      = list_sColumn[15]
        cClipData.sAssembly     = list_sColumn[40]

        # Filter out BAM or other file formats
        if cClipData.sFormat   != 'bed narrowPeak': continue
        if cClipData.sAssembly != 'hg19':           continue

        list_sOutput.append(cClipData)
    # loop END: i, sReadLine
    InFile.close()

    # VS Check
    if not list_sOutput:
        sys.exit('Invalid List : eCLIP_parse_eclip_meta : list_sOutput size= %d' % (len(list_sOutput)))

    return list_sOutput
## endregion

## region class cPeakData
class cPeakData:
    def __init__(self):
        self.sChrID         = ''
        self.nStartPos      = 0
        self.nEndPos        = 0
        self.nCenterPos     = 0
        self.nPeakSize      = 0
        self.sSample        = ''
        self.nScore         = 0
        self.sStrand        = ''
        self.sLocale        = ''  #ORF, 5UTR, 3UTR, or Intergenic
        self.list_cVCF      = []

        ## Extra
        self.sGeneSym           = ''
        self.sNMID              = ''
        self.list_fDepth_CLIP   = []
        self.list_fDepth_mRNA   = []
        self.list_nPeakPos      = []
        self.nWigStart          = 0
        self.nWigEnd            = 0
        # def END: __init__

def cPeakData_parse_bedfile(sInFile, dict_cRefSeq, bBuff = 0):

    dict_cRef      = {}
    for sNMID in dict_cRefSeq:
        sKey = '%s,%s' % (dict_cRefSeq[sNMID].sChrID, dict_cRefSeq[sNMID].sStrand)
        if sKey not in dict_cRef:
            dict_cRef[sKey] = []
        dict_cRef[sKey].append(dict_cRefSeq[sNMID])
    #loop END: sNMID

    list_sOutput   = []
    InFile         = open(sInFile, 'r')
    for sReadLine in InFile:
        # File Format:
        # Column:        0       | 1         | 2        | 3               | 4      | 5       |
        # ID:            sChrID  | nStartPos | nEndPos  | sSample         | nScore | sStrand |
        # Example:       chr14   | 56879239  | 56879435 | ILF3_K562_rep02 | 1000   | -       |
        # Column:        6       | 7         | 8        | 9
        # ID:            fSignal | fPvalue   | Not used | Not used
        # Example:       1.29065 | 0.198802  | -1       | -1

        list_sColumn     = sReadLine.strip('\n').split('\t')

        cPeak            = cPeakData()
        cPeak.sChrID     = list_sColumn[0]
        cPeak.nStartPos  = int(list_sColumn[1])
        cPeak.nEndPos    = int(list_sColumn[2])

        assert cPeak.nStartPos < cPeak.nEndPos

        cPeak.nPeakSize  = cPeak.nEndPos - cPeak.nStartPos
        cPeak.nCenterPos = cPeak.nEndPos - int(cPeak.nPeakSize / 2)

        cPeak.sSample    = list_sColumn[3]
        cPeak.nScore     = int(list_sColumn[4])
        cPeak.sStrand    = list_sColumn[5]
        cPeak.fSignal    = list_sColumn[6]
        cPeak.fPvalue    = list_sColumn[7]
        cPeak.sNotUsed1  = list_sColumn[8]
        cPeak.sNotUsed2  = list_sColumn[9]

        # Determine Peak Genomic Location
        if not bBuff: determine_peak_locale(cPeak, dict_cRef)
        list_sOutput.append(cPeak)
    # loop END: sReadLine
    InFile.close()

    #VS Check:
    if not list_sOutput:
        sys.exit('ERROR: Invalid list_sOutput Size= %d' % (len(list_sOutput)))

    return list_sOutput
#def END: cPeakData_parse_bedfile

def cPeakData_parse_bedfile_alt(sInFile, dict_cRefSeq):

    dict_cRef      = {}
    for sNMID in dict_cRefSeq:
        sKey = '%s,%s' % (dict_cRefSeq[sNMID].sChrID, dict_cRefSeq[sNMID].sStrand)
        if sKey not in dict_cRef:
            dict_cRef[sKey] = []
        dict_cRef[sKey].append(dict_cRefSeq[sNMID])
    #loop END: sNMID

    list_sOutput   = []
    InFile         = open(sInFile, 'r')
    for sReadLine in InFile:
        # File Format:
        # Column:        0       | 1         | 2        | 3        | 4      | 5       |
        # ID:            sChrID  | nStartPos | nEndPos  | sSample  | nScore | sStrand |
        # Example:       chr14   | 56879239  | 56879435 | Gzmb     | 1000   | -       |
        # Column:        6       | 7         |
        # ID:            fSignal | fPvalue   |
        # Example:       1.29065 | 0.198802  |

        list_sColumn     = sReadLine.strip('\n').split('\t')

        cPeak            = cPeakData()
        cPeak.sChrID     = list_sColumn[0]
        cPeak.nStartPos  = int(list_sColumn[1])
        cPeak.nEndPos    = int(list_sColumn[2])

        assert cPeak.nStartPos < cPeak.nEndPos

        cPeak.nPeakSize  = cPeak.nEndPos - cPeak.nStartPos
        cPeak.nCenterPos = cPeak.nEndPos - int(cPeak.nPeakSize / 2)

        cPeak.list_sSample    = list_sColumn[3].split(',')


        cPeak.list_nScore     = list_sColumn[4].split(',')
        cPeak.sStrand    = list_sColumn[5].split(',')[0]
        cPeak.list_sStrand    = list_sColumn[5].split(',')[0]
        cPeak.list_fSignal    = list_sColumn[6].split(',')
        cPeak.list_fPvalue    = list_sColumn[7].split(',')

        determine_peak_locale(cPeak, dict_cRef)

        list_sOutput.append(cPeak)
    # loop END: sReadLine
    InFile.close()

    #VS Check:
    if not list_sOutput:
        sys.exit('ERROR: Invalid list_sOutput Size= %d' % (len(list_sOutput)))

    return list_sOutput
#def END: cPeakData_parse_bedfile_alt


def determine_peak_locale (cPeak, dict_cRef):

    if cPeak.sChrID == 'chrM': cPeak.sLocale = 'mitochondrial'

    else:
        list_cRef = dict_cRef['%s,%s' % (cPeak.sChrID, cPeak.sStrand)]
        list_cRef = [cRef for cRef in list_cRef if cRef.nTxStart <= cPeak.nCenterPos <= cRef.nTxEnd]

        if len(list_cRef) == 0:  #Peak CenterPos not within transcript
            cPeak.sGeneSym = 'intergenic'
            cPeak.sNMID    = 'intergenic'
            cPeak.sLocale  = 'intergenic'

        elif len(list_cRef) == 1:
            cRef            = list_cRef[0]
            cPeak.sGeneSym  = cRef.sGeneSym
            cPeak.sNMID     = cRef.sNMID

            if cPeak.nCenterPos < cRef.nCdsStart:
                cPeak.sLocale = '5UTR' if cPeak.sStrand == '+' else '3UTR'

            elif cPeak.nCenterPos > cRef.nCdsEnd:
                cPeak.sLocale = '3UTR' if cPeak.sStrand == '+' else '5UTR'

            else:
                for nStart, nEnd in zip(cRef.list_nORFS, cRef.list_nORFE):

                    if nStart <= cPeak.nCenterPos <= nEnd:
                        cPeak.sLocale = 'exonic'
                        break
                    else: cPeak.sLocale = 'intronic'
                #loop END: nStart, nEnd
            #if END:
        elif len(list_cRef) == 2: # More than 1 gene (nested gene)

            # Prioritize by UTRs > ORF > Intron
            for i, cRef in enumerate(list_cRef):

                if cPeak.nCenterPos < cRef.nCdsStart:
                    cPeak.sGeneSym = cRef.sGeneSym
                    cPeak.sNMID    = cRef.sNMID
                    cPeak.sLocale  = '5UTR' if cPeak.sStrand == '+' else '3UTR'
                    break
                elif cPeak.nCenterPos > cRef.nCdsEnd:
                    cPeak.sGeneSym = cRef.sGeneSym
                    cPeak.sNMID    = cRef.sNMID
                    cPeak.sLocale  = '3UTR' if cPeak.sStrand == '+' else '5UTR'
                    break

                else:
                    for nStart, nEnd in zip(cRef.list_nORFS, cRef.list_nORFE):

                        if nStart <= cPeak.nCenterPos <= nEnd:
                            cPeak.sGeneSym = cRef.sGeneSym
                            cPeak.sNMID    = cRef.sNMID
                            cPeak.sLocale  = 'exonic'
                            break
                        else:
                            if i == 0:
                                cPeak.sGeneSym = cRef.sGeneSym
                                cPeak.sNMID    = cRef.sNMID
                                cPeak.sLocale  = 'intronic'

                    #loop END: nStart, nEnd
                #if END:
            # loop END: cRef
        #if END:
    #if END:
#def END: determine_peak_local
## endregion

## region class cRefSeq

class cRefSeq:
    def __init__(self):
        # Initalization of instances variables
        self.sGeneSym       = ''
        self.sNMID          = ''
        self.nChrID         = 0
        self.sChrID         = 0
        self.sStrand        = ''
        self.nTxStart       = 0
        self.nTxEnd         = 0
        self.nCdsStart      = 0
        self.nCdsEnd        = 0
        self.nExonCount     = 0

        self.sORFSeq        = ''
        self.s5UTRSeq       = ''
        self.s3UTRSeq       = 'NULL'

        self.sORFStart      = 0
        self.sORFEnd        = 0
        self.s5UTRStart     = 0
        self.s5UTREnd       = 0
        self.s3UTRStart     = 0
        self.s3UTREnd       = 0

        self.n5UTRLen       = 0
        self.n3UTRLen       = 0
        self.nORFLen        = 0

        # Exon Start and End are lists with position nExonStart[0] with nExonEnd[0[ and so fourth
        self.list_nExonS    = []
        self.list_nExonE    = []
        self.list_n5UTRS    = []
        self.list_n5UTRE    = []
        self.list_nORFS     = []
        self.list_nORFE     = []
        self.list_n3UTRS    = []
        self.list_n3UTRE    = []
        self.list_sCapture  = []

    # def END: __init__

def cRef_parse_refflat_line(sInFile, sOutFile):

    if not os.path.isfile(sOutFile):
        print('Pickling Refflat', sOutFile)
        cGenome         = cFasta(dict_sGENOME[sGENOME])
        dict_sOutput    = {}
        InFile          = open(sInFile, 'r')
        nCnt            = 0
        nCnt1           = 0
        for sReadLine in InFile:

            list_sColumn = sReadLine.strip('\n').split('\t')
            cRef                = cRefSeq()
            cRef.sGeneSym       = list_sColumn[0].upper()
            cRef.sNMID          = list_sColumn[1]

            #if cRef.sNMID != 'NM_016200': continue
            #if 'ALDH3B1' not in cRef.sGeneSym: continue

            cRef.sChrID         = list_sColumn[2]
            cRef.nChrID         = process_chromo_ID(list_sColumn[2])
            cRef.sStrand        = list_sColumn[3]
            cRef.nTxStart       = int(list_sColumn[4])
            cRef.nTxEnd         = int(list_sColumn[5])
            cRef.nCdsStart      = int(list_sColumn[6])
            cRef.nCdsEnd        = int(list_sColumn[7])
            cRef.nExonCount     = int(list_sColumn[8])
            cRef.list_nExonS    = [int(sStartPos) for sStartPos in list_sColumn[9][:-1].split(',')]
            cRef.list_nExonE    = [int(nEndPos)   for nEndPos   in list_sColumn[10][:-1].split(',')]
            assert len(cRef.list_nExonS) == len(cRef.list_nExonE)
            #print(cRef.sGeneSym, cRef.sStrand, cRef.nCdsStart, cRef.nCdsEnd)
            cRef_determine_seq (cRef, cGenome)

            cRef.nIntronSize =  (cRef.nTxEnd - cRef.nTxStart) - (len(cRef.s5UTRSeq) + len(cRef.sORFSeq) + len(cRef.s3UTRSeq))

            if cRef.sGeneSym not in dict_sOutput:
                dict_sOutput[cRef.sGeneSym] = ''
            dict_sOutput[cRef.sGeneSym] = cRef

        #loop END: sReadLine
        InFile.close()
        #VS Check
        if not dict_sOutput:
            sys.exit('Invalid Dictionary : cRef_parse_refflat_line : dict_sOutput size= %d' % len(dict_sOutput))

        #Pickle dump
        OutFile = open(sOutFile, 'wb')
        pickle.dump(dict_sOutput, OutFile)
        OutFile.close()

        return dict_sOutput

    else:
        print('Opening Pickled', sOutFile)
        OutFile      = open(sOutFile, 'rb')
        dict_sOutput =  pickle.load(OutFile)
        OutFile.close()
        # VS Check
        if not dict_sOutput:
            sys.exit('Invalid Dictionary : cRef_parse_refflat_line : dict_sOutput size= %d' % len(dict_sOutput))
        return dict_sOutput
    #if END:
#def END: cRef_parse_refflat_line

def cRef_parse_refflat_line_v2(sInFile, list_sCapture, bCapture = 0):

    sOutDir         = '%s/pickles' % sBASE_DIR
    os.makedirs(sOutDir, exist_ok=True)

    dict_sCapGene   = list_sCapture[0]
    dict_sCapID     = list_sCapture[1]

    if bCapture: sRefSeqFile = '%s/refseq.capture3.dat' % sOutDir
    else:        sRefSeqFile = '%s/refseq.dat' % sOutDir

    if not os.path.isfile(sRefSeqFile):
        print('Pickling Refflat', sRefSeqFile)
        cGenome         = cFasta(dict_sGENOME[sGENOME])
        dict_sOutput    = {}
        InFile          = open(sInFile, 'r')
        nCnt            = 0
        nCnt1           = 0
        for sReadLine in InFile:

            list_sColumn = sReadLine.strip('\n').split('\t')
            cRef                = cRefSeq()
            cRef.sGeneSym       = list_sColumn[0].upper()
            cRef.sNMID          = list_sColumn[1]

            #if cRef.sNMID != 'NM_016200': continue
            #if 'ALDH3B1' not in cRef.sGeneSym: continue

            cRef.sChrID         = list_sColumn[2]
            cRef.nChrID         = process_chromo_ID(list_sColumn[2])
            cRef.sStrand        = list_sColumn[3]
            cRef.nTxStart       = int(list_sColumn[4])
            cRef.nTxEnd         = int(list_sColumn[5])
            cRef.nCdsStart      = int(list_sColumn[6])
            cRef.nCdsEnd        = int(list_sColumn[7])
            cRef.nExonCount     = int(list_sColumn[8])
            cRef.list_nExonS    = [int(sStartPos) for sStartPos in list_sColumn[9][:-1].split(',')]
            cRef.list_nExonE    = [int(nEndPos)   for nEndPos   in list_sColumn[10][:-1].split(',')]
            assert len(cRef.list_nExonS) == len(cRef.list_nExonE)
            #print(cRef.sGeneSym, cRef.sStrand, cRef.nCdsStart, cRef.nCdsEnd)

            dict_sPos           = {}
            if bCapture:
                # GeneSym and NMID do not match on same chrID
                try: dict_sPos.update(dict_sCapGene[cRef.sGeneSym])
                except KeyError: dict_sPos.update({})

                try: dict_sPos.update(dict_sCapID[cRef.sNMID])
                except KeyError: dict_sPos.update({})

                #print(dict_sPos)
                #sys.exit()

                if not dict_sPos:
                    nCnt += 1
                    continue # 246/18625 genes not found in capture file

                bCovered = cRef_determine_capture_coverage (cRef, dict_sPos)
                if not bCovered:
                    nCnt1 += 1
                    continue
            #if END: bCapture

            cRef_determine_seq (cRef, cGenome, bCapture)

            if cRef.sNMID not in dict_sOutput:
                dict_sOutput[cRef.sNMID] = ''
            dict_sOutput[cRef.sNMID] = cRef

            nCnt += sum([e-s +1 for s ,e in zip(cRef.list_n5UTRS, cRef.list_n5UTRE)])
            nCnt1 += sum([e-s +1 for s ,e in zip(cRef.list_n3UTRS, cRef.list_n3UTRE)])



        #loop END: sReadLine
        InFile.close()
        #VS Check
        if not dict_sOutput:
            sys.exit('Invalid Dictionary : cRef_parse_refflat_line : dict_sOutput size= %d' % len(dict_sOutput))

        #Pickle dump
        OutFile = open(sRefSeqFile, 'wb')
        pickle.dump(dict_sOutput, OutFile)
        OutFile.close()

        print('Dumped', len(dict_sOutput), 'NotFound', nCnt, 'NoCover', nCnt1)

        return dict_sOutput

    else:
        print('Opening Pickled', sRefSeqFile)
        dict_sOutput =  pickle.load(open(sRefSeqFile, 'rb'))

        # VS Check
        if not dict_sOutput:
            sys.exit('Invalid Dictionary : cRef_parse_refflat_line : dict_sOutput size= %d' % len(dict_sOutput))
        return dict_sOutput
    #if END:
#def END: cRef_parse_refflat_line

def cRef_determine_capture_coverage (cRef, dict_sPos):

    # Update exon positions with coverage overlap
    list_nCoverage   = []
    list_n5UTR       = []
    list_n3UTR       = []
    list_nIntrons    = []
    list_nIntergenic = []

    list_sCapture    = dict_sPos[cRef.sChrID] #VS-Check
    list_sCapture    = [[int(sPosKey.split(',')[0]), int(sPosKey.split(',')[1])]
                        for sPosKey in list_sCapture]
    list_sCapture    = sorted(list_sCapture, key=lambda e:e[0])

    for i, (nExonS, nExonE) in enumerate(zip(cRef.list_nExonS, cRef.list_nExonE)):
        #print(red('Exons', 1), i, nExonS, nExonE)

        for nCapStart, nCapEnd in list_sCapture:

            if nCapEnd < nExonS:   continue
            if nCapStart > nExonE: continue
            nNewStart = max(nCapStart, nExonS)
            nNewEnd   = min(nCapEnd, nExonE)
            #print(green('Capture',1), nCapStart, nCapEnd)
            #print(blue('Overlap',1), nNewStart, nNewEnd)

            if nCapStart < cRef.nCdsStart: # 5UTR or if long enough intergenic in the 5' region

                if nNewEnd > max(nNewStart, cRef.nCdsStart):
                    list_nCoverage.append([max(nNewStart, cRef.nCdsStart), nNewEnd])

                if nNewStart < cRef.nCdsStart:
                    list_n5UTR.append([nNewStart, min(nNewEnd, cRef.nCdsStart)])

                if nCapStart < cRef.nTxStart:
                    list_nIntergenic.append([nCapStart, cRef.nTxStart])

                if nCapEnd > nExonE:
                    list_nIntrons.append([nExonE, nCapEnd])

            elif nCapEnd > cRef.nCdsEnd: # 3UTR or if long enough intergenic in the 3' region

                if i != 0:
                    if nNewStart < min(nNewEnd, cRef.nCdsEnd):
                        list_nCoverage.append([nNewStart, min(nNewEnd, cRef.nCdsEnd)])

                if nNewEnd > cRef.nCdsEnd:
                    list_n3UTR.append([max(cRef.nCdsEnd, nCapStart), nNewEnd])

                if nCapEnd > cRef.nTxEnd:
                    list_nIntergenic.append([cRef.nTxEnd, nCapEnd])

                if nCapStart < nExonS:
                    list_nIntrons.append([nCapStart, nExonS])

            else:
                list_nCoverage.append([nNewStart, nNewEnd])

                if nCapStart < nExonS:
                    list_nIntrons.append([nCapStart, nExonS])
                if nCapEnd > nExonE:
                    list_nIntrons.append([nExonE, nCapEnd])
            #if END:
        # loop END: sPosKey
    #loop END: i, (nExonS, nExonE)
    '''
    print(list_nIntergenic)
    for s,e in list_nIntergenic:
        print(s,e, e-s)

    nCaptureStart    = list_sCapture[0][0]
    nCaptureEnd      = list_sCapture[-1][1]

    nInter5          = cRef.nTxStart - nCaptureStart
    nInter3          = nCaptureEnd - cRef.nTxEnd
    if nInter5 > 0 and nInter3 > 0:
        print(cRef.nTxStart, cRef.nTxEnd)
        print(nCaptureStart, nCaptureEnd)
        print(nInter5, nInter3)
        sys.exit()
    list_nCoverage  = sorted(list_nCoverage, key=lambda e:e[0])

    print('list_nIntergenic', list_nIntergenic)
    print('list_nIntrons', list_nIntrons)
    print('list_n5UTR', list_n5UTR)
    print('list_n3UTR', list_n3UTR)
    print('list_nCoverage', len(list_nCoverage))

    for e in list_n5UTR:print('5', e, e[1]-e[0])
    for e in list_nCoverage:print('e', e, e[1]-e[0])
    for e in list_n3UTR:print('3', e, e[1]-e[0])
    for e in list_nIntrons:print('i', e, e[1]-e[0])
    sys.exit()
    '''

    if not list_nCoverage: # No overlap
        return 0
    else:
        cRef.list_nIntergenic = list_nIntergenic
        cRef.list_nIntrons    = list_nIntrons
        cRef.list_n5UTRS      = [nStart for nStart, nEnd in list_n5UTR if nStart != nEnd]
        cRef.list_n5UTRE      = [nEnd   for nStart, nEnd in list_n5UTR if nStart != nEnd]
        cRef.list_n3UTRS      = [nStart for nStart, nEnd in list_n3UTR if nStart != nEnd]
        cRef.list_n3UTRE      = [nEnd   for nStart, nEnd in list_n3UTR if nStart != nEnd]
        cRef.list_nExonS      = [nStart for nStart, nEnd in list_nCoverage if nStart != nEnd]
        cRef.list_nExonE      = [nEnd   for nStart, nEnd in list_nCoverage if nStart != nEnd]
        cRef.nCaptureS        = cRef.list_nExonS[0]
        cRef.nCaptureE        = cRef.list_nExonE[-1]
        if cRef.nCaptureS > cRef.nCdsEnd: return 0
        elif cRef.nCaptureE < cRef.nCdsStart: return 0
        else: return 1
#def END: cRef_check_capture_coverage

def cRef_determine_seq (cRef, cGenome):

    nCdsStart     = cRef.nCdsStart
    nCdsEnd       = cRef.nCdsEnd
    s5UTRSeq      = ''
    sORFSeq       = ''
    s3UTRSeq      = ''

    list_n5UTRS   = []
    list_n5UTRE   = [nCdsStart]
    list_nORFS    = [nCdsStart]
    list_nORFE    = [nCdsEnd]
    list_n3UTRS   = [nCdsEnd]
    list_n3UTRE   = []

    # Divide up the exon start positions into three lists 5', ORF, and 3'
    for nStartPos in cRef.list_nExonS:

        if nStartPos <= nCdsStart:
            list_n5UTRS.append(nStartPos) # Positions before CdsStart

        elif nStartPos >= nCdsStart and nStartPos < nCdsEnd:
            list_nORFS.append(nStartPos) # Positions Between CdsStart and CdsEnd
            list_nORFS = sorted(list_nORFS)

        else:
            list_n3UTRS.append(nStartPos)
            list_n3UTRS = sorted(list_n3UTRS) # Positions after CdsEnd
        #if END:
    #loop END: nStartPos

    # Divide up the exon end positions into three lists 5', ORF, and 3'
    for nEndPos in cRef.list_nExonE:

        if nEndPos <= nCdsStart:
            list_n5UTRE.append(nEndPos)
            list_n5UTRE = sorted(list_n5UTRE)

        elif nEndPos >= nCdsStart and nEndPos < nCdsEnd:
            list_nORFE.append(nEndPos)
            list_nORFE = sorted(list_nORFE)

        else:
            list_n3UTRE.append(nEndPos)
        #if END:
    #loop END: nEndPos

    for nStartPos, nEndPos in zip(list_n5UTRS, list_n5UTRE):
        if nStartPos == nEndPos: continue
        #print('5UTR', nStartPos, nEndPos, nEndPos-nStartPos)
        s5UTRSeq += cGenome.fetch(cRef.sChrID, cRef.sGeneSym, nStartPos, nEndPos, cRef.sStrand)

    for nStartPos, nEndPos in zip(list_nORFS, list_nORFE):
        #print('ORF', nStartPos, nEndPos, nEndPos-nStartPos)
        sORFSeq += cGenome.fetch(cRef.sChrID, cRef.sGeneSym, nStartPos, nEndPos, cRef.sStrand)

    for nStartPos, nEndPos in zip(list_n3UTRS, list_n3UTRE):
        if nStartPos == nEndPos: continue
        #print('3UTR', nStartPos, nEndPos, nEndPos-nStartPos)
        s3UTRSeq += cGenome.fetch(cRef.sChrID, cRef.sGeneSym, nStartPos, nEndPos, cRef.sStrand)

    # Compile each segement by range slicing the original chromosome sequence
    if cRef.sStrand == '-':
        # For '-' strand, switch the 5' and the 3'
        cRef.list_n5UTRS = list_n3UTRS
        cRef.list_n5UTRE = list_n3UTRE
        cRef.list_n3UTRS = list_n5UTRS
        cRef.list_n3UTRE = list_n5UTRE
        cRef.list_nORFS  = list_nORFS
        cRef.list_nORFE  = list_nORFE
        cRef.s5UTRSeq    = reverse_complement(s3UTRSeq)
        cRef.s3UTRSeq    = reverse_complement(s5UTRSeq)
        cRef.sORFSeq     = reverse_complement(sORFSeq)
    else:
        cRef.list_n5UTRS = list_n5UTRS
        cRef.list_n5UTRE = list_n5UTRE
        cRef.list_n3UTRS = list_n3UTRS
        cRef.list_n3UTRE = list_n3UTRE
        cRef.list_nORFS  = list_nORFS
        cRef.list_nORFE  = list_nORFE
        cRef.s5UTRSeq    = s5UTRSeq
        cRef.s3UTRSeq    = s3UTRSeq
        cRef.sORFSeq     = sORFSeq
    #if END:
    #print('5UTR', len(cRef.s5UTRSeq))
    #print('ORF',  len(cRef.sORFSeq))
    #print('3UTR', len(cRef.s3UTRSeq))
    #sys.exit()
#def END: cRef_determine_seq

def cRef_update_seqinfo (cRef, cGenome):

    s5UTRSeq      = ''
    sORFSeq       = ''
    s3UTRSeq      = ''

    list_n5UTRS   = cRef.list_n5UTRS
    list_n5UTRE   = cRef.list_n5UTRE
    list_nORFS    = cRef.list_nExonS
    list_nORFE    = cRef.list_nExonE
    list_n3UTRS   = cRef.list_n3UTRS
    list_n3UTRE   = cRef.list_n3UTRE


    for nStartPos, nEndPos in zip(list_n5UTRS, list_n5UTRE):
        if nStartPos == nEndPos: continue
        #print('5UTR', nStartPos, nEndPos, nEndPos-nStartPos)
        s5UTRSeq += cGenome.fetch(cRef.sChrID, cRef.sGeneSym, nStartPos, nEndPos, cRef.sStrand)

    for nStartPos, nEndPos in zip(list_nORFS, list_nORFE):
        #print('ORF', nStartPos, nEndPos, nEndPos-nStartPos)
        sORFSeq += cGenome.fetch(cRef.sChrID, cRef.sGeneSym, nStartPos, nEndPos, cRef.sStrand)

    for nStartPos, nEndPos in zip(list_n3UTRS, list_n3UTRE):
        if nStartPos == nEndPos: continue
        #print('3UTR', nStartPos, nEndPos, nEndPos-nStartPos)
        s3UTRSeq += cGenome.fetch(cRef.sChrID, cRef.sGeneSym, nStartPos, nEndPos, cRef.sStrand)

    # Compile each segement by range slicing the original chromosome sequence
    if cRef.sStrand == '-':
        # For '-' strand, switch the 5' and the 3'
        cRef.list_n5UTRS = list_n3UTRS
        cRef.list_n5UTRE = list_n3UTRE
        cRef.list_n3UTRS = list_n5UTRS
        cRef.list_n3UTRE = list_n5UTRE
        cRef.list_nORFS  = list_nORFS
        cRef.list_nORFE  = list_nORFE
        cRef.s5UTRSeq    = reverse_complement(s3UTRSeq)
        cRef.s3UTRSeq    = reverse_complement(s5UTRSeq)
        cRef.sORFSeq     = reverse_complement(sORFSeq)
    else:
        cRef.list_n5UTRS = list_n5UTRS
        cRef.list_n5UTRE = list_n5UTRE
        cRef.list_n3UTRS = list_n3UTRS
        cRef.list_n3UTRE = list_n3UTRE
        cRef.list_nORFS  = list_nORFS
        cRef.list_nORFE  = list_nORFE
        cRef.s5UTRSeq    = s5UTRSeq
        cRef.s3UTRSeq    = s3UTRSeq
        cRef.sORFSeq     = sORFSeq
    #if END:
    #print('5UTR', len(cRef.s5UTRSeq))
    #print('ORF',  len(cRef.sORFSeq))
    #print('3UTR', len(cRef.s3UTRSeq))
    #sys.exit()
#def END: cRef_determine_seq

def cRef_determine_seq_mod (cRef, cGenome, bCover):

    nCdsStart     = cRef.nCdsStart if cRef.nCdsStart > cRef.nCoverStart and cRef.nCdsStart < cRef.nCoverEnd else cRef.nCoverStart
    nCdsEnd       = cRef.nCdsEnd if cRef.nCdsEnd < cRef.nCoverEnd and cRef.nCdsEnd > cRef.nCoverStart else cRef.nCoverEnd
    s5UTRSeq      = ''
    sORFSeq       = ''
    s3UTRSeq      = ''

    '''
    print('CDSs', nCdsStart, nCdsEnd, nCdsEnd - nCdsStart)
    print('bCover', bCover)
    for i, (nExonS, nExonE) in enumerate(zip(cRef.list_nExonS, cRef.list_nExonE)):
        print('NewExon', i+1, nExonS, nExonE, nExonE-nExonS)
    '''

    if bCover == 1: # Coverage only in 3'UTR
        list_n5UTRS   = []
        list_n5UTRE   = []
        list_nORFS    = []
        list_nORFE    = []
        list_n3UTRS   = cRef.list_nExonS
        list_n3UTRE   = cRef.list_nExonE
    elif bCover == 2: # Coverage only in 5'UTR
        list_n5UTRS   = cRef.list_nExonS
        list_n5UTRE   = cRef.list_nExonE
        list_nORFS    = []
        list_nORFE    = []
        list_n3UTRS   = []
        list_n3UTRE   = []
    else:
        list_n5UTRS   = []
        list_n5UTRE   = [nCdsStart]
        list_nORFS    = [nCdsStart]
        list_nORFE    = [nCdsEnd]
        list_n3UTRS   = [nCdsEnd]
        list_n3UTRE   = []
        # Divide up the exon start positions into three lists 5', ORF, and 3'
        for nStartPos in cRef.list_nExonS:

            if nStartPos <= nCdsStart:
                list_n5UTRS.append(nStartPos) # Positions before CdsStart

            elif nStartPos >= nCdsStart and nStartPos < nCdsEnd:
                list_nORFS.append(nStartPos) # Positions Between CdsStart and CdsEnd
                list_nORFS = sorted(list_nORFS)

            else:
                list_n3UTRS.append(nStartPos)
                list_n3UTRS = sorted(list_n3UTRS) # Positions after CdsEnd
            #if END:
        #loop END: nStartPos

        # Divide up the exon end positions into three lists 5', ORF, and 3'
        for nEndPos in cRef.list_nExonE:

            if nEndPos <= nCdsStart:
                list_n5UTRE.append(nEndPos)
                list_n5UTRE = sorted(list_n5UTRE)

            elif nEndPos >= nCdsStart and nEndPos < nCdsEnd:
                list_nORFE.append(nEndPos)
                list_nORFE = sorted(list_nORFE)

            else:
                list_n3UTRE.append(nEndPos)
            #if END:
        #loop END: nEndPos
    #if END:

    for nStartPos, nEndPos in zip(list_n5UTRS, list_n5UTRE):
        if nStartPos == nEndPos: continue
        #print('5UTR', nStartPos, nEndPos, nEndPos-nStartPos)
        s5UTRSeq += cGenome.fetch(cRef.sChrID, cRef.sGeneSym, nStartPos, nEndPos, cRef.sStrand)

    for nStartPos, nEndPos in zip(list_nORFS, list_nORFE):
        #print('ORF', nStartPos, nEndPos, nEndPos-nStartPos)
        sORFSeq += cGenome.fetch(cRef.sChrID, cRef.sGeneSym, nStartPos, nEndPos, cRef.sStrand)

    for nStartPos, nEndPos in zip(list_n3UTRS, list_n3UTRE):
        if nStartPos == nEndPos: continue
        #print('3UTR', nStartPos, nEndPos, nEndPos-nStartPos)
        s3UTRSeq += cGenome.fetch(cRef.sChrID, cRef.sGeneSym, nStartPos, nEndPos, cRef.sStrand)

    # Compile each segement by range slicing the original chromosome sequence
    if cRef.sStrand == '-':
        # For '-' strand, switch the 5' and the 3'
        cRef.list_n5UTRS = list_n3UTRS
        cRef.list_n5UTRE = list_n3UTRE
        cRef.list_n3UTRS = list_n5UTRS
        cRef.list_n3UTRE = list_n5UTRE
        cRef.list_nORFS  = list_nORFS
        cRef.list_nORFE  = list_nORFE
        cRef.s5UTRSeq    = reverse_complement(s3UTRSeq)
        cRef.s3UTRSeq    = reverse_complement(s5UTRSeq)
        cRef.sORFSeq     = reverse_complement(sORFSeq)
    else:
        cRef.list_n5UTRS = list_n5UTRS
        cRef.list_n5UTRE = list_n5UTRE
        cRef.list_n3UTRS = list_n3UTRS
        cRef.list_n3UTRE = list_n3UTRE
        cRef.list_nORFS  = list_nORFS
        cRef.list_nORFE  = list_nORFE
        cRef.s5UTRSeq    = s5UTRSeq
        cRef.s3UTRSeq    = s3UTRSeq
        cRef.sORFSeq     = sORFSeq
    #if END:
    #print('5UTR', len(cRef.s5UTRSeq))
    #print('ORF',  len(cRef.sORFSeq))
    #print('3UTR', len(cRef.s3UTRSeq))
#def END: cRef_determine_seq

## region class c3UTRs

class c3UTRs:
    def __init__(self):
        self.sStrand    = 'NULL'
        self.sMAID      = 'NULL'
        self.sMirSeq    = 'NULL'
        self.sNMID      = 'NULL'
        self.sGeneSym   = 'NULL'
        self.n3UTRLen   = 0
        self.s3UTRSeq   = 'NULL'

        #For Compare Repression
        self.bSTFlag     = 'NULL'
        self.fFC         = 0.0
        self.fExpr       = 0.0
        self.dict_nSTCnt = {}
        self.dict_nSTPos = {}

        #For Multi-Linear Regression
        self.fTA         = {}
        self.fSPS        = {}
        self.fLocalAU    = {}

        # From RefSeq
        self.list_n5UTRS    = []
        self.list_n5UTRE    = []
        self.list_nORFS     = []
        self.list_nORFE     = []
        self.list_n3UTRS    = []
        self.list_n3UTRE    = []
    #def END: __init__

def c3UTRs_load_3UTRS (sInFileDir,  b3UTR_SIZE_FILTER=False ):

    dict_c3UTRs  = {}

    InFile       = open(sInFileDir, 'r')
    for sReadLine in InFile:
        # File Format
        # Column Number:     | 0                     | 1        | 2           |
        # Column Description:| sStrand               | n3UTRLen | s3UTRSeq    |
        # Column Example:    | -|NM_001003806|TARP|* | 455      | CAGACGGT... |

        list_sColumn = sReadLine.strip('\n').split('\t')

        # V-S Check:
        if len(list_sColumn) != n3UTR_FILE_COLUMN_SIZE:
            sys.exit('ERROR: 3UTR Text File Column Size= %d\tContent= %s' % (len(list_sColumn), list_sColumn))

        c3UTR            = c3UTRs()
        c3UTR.sStrand    =     list_sColumn[0].split('|')[0]
        c3UTR.sNMID      =     list_sColumn[0].split('|')[1]
        c3UTR.sGeneSym   =     list_sColumn[0].split('|')[2].upper()
        c3UTR.n3UTRLen   = int(list_sColumn[1])
        c3UTR.s3UTRSeq   =     list_sColumn[2].upper().replace('T', 'U')

        if c3UTR.n3UTRLen > 0:
            c3UTR.fAUScore   = (c3UTR.s3UTRSeq.count('A') + c3UTR.s3UTRSeq.count('U') + c3UTR.s3UTRSeq.count('T')) / c3UTR.n3UTRLen
        else:
            c3UTR.fAUScore   = 0.5


        if b3UTR_SIZE_FILTER == True:
            if c3UTR.n3UTRLen > nMAX_3UTR_SIZE: continue

        if c3UTR.sNMID not in dict_c3UTRs:
            dict_c3UTRs[c3UTR.sNMID] = ''

        dict_c3UTRs[c3UTR.sNMID] = c3UTR
        #list_cCandids.append(cCandid)
    #loop END: sReadLine

    InFile.close()

    #V-S Check: Dictionary Size
    if not dict_c3UTRs:
        sys,exit('ERROR: dict_c3UTRs Size= %d' % len(dict_c3UTRs))


    return dict_c3UTRs
#def END: c3UTRs_load_3UTRS

def c3UTRs_normalize_FC (c3UTR, fAvgFC):
    #c3UTR     = copy.deepcopy(c3UTR)
    c3UTR.fFC = c3UTR.fFC - fAvgFC
    return c3UTR
#def END: c3UTRs_normalize_FC

def c3UTRs_assign_TA_and_SPS (c3UTR, dict_TA_SPS):

    s7m8Seq               = determine_site_types(c3UTR.sMirSeq)['7mer-m8']
    c3UTR.fTA, c3UTR.fSPS = dict_TA_SPS[s7m8Seq]

    return c3UTR
#def END: c3UTRs_assign_TA_and_SPS

## endregion


## endregion

## region class cVCFData

class cVCFData:
    def __init__(self):
        self.sPatID         = ''
        self.sChrID         = ''
        self.nPos           = 0
        self.sDBSNP_ID      = ''
        self.sRefNuc        = ''
        self.sAltNuc        = ''
        self.fQual          = 0.0
        self.sFilter        = ''
        self.sInfo          = ''
        self.sFormat        = ''
        self.list_sMisc     = []

        #Extra
        self.nClusterID     = 0

    #def END: __int__


def cVCF_parse_vcf_files (sVCFFile, bClustered=False):
    if not os.path.isfile(sVCFFile):
        sys.exit('File Not Found %s' % sVCFFile)

    list_sOutput = []
    InFile       = open(sVCFFile, 'r')

    for sReadLine in InFile:
        # File Format
        # Column Number:     | 0       | 1        | 2          | 3       | 4
        # Column Description:| sChrID  | nPos     | sDBSNP_ID  | sRefNuc | sAltNuc
        # Column Example:    | chr13   | 32906558 | rs79483201 | T       | A
        # Column Number:     | 5       | 6        | 7          | 8              | 9./..
        # Column Description:| fQual   | sFilter  | sInfo      | sFormat        | sSampleIDs
        # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to sFormat
        ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph

        if sReadLine.startswith('#'): continue  # SKIP Information Headers
        list_sColumn        = sReadLine.strip('\n').split('\t')

        cVCF                = cVCFData()
        if sSEQDATA == 'WXS':
            cVCF.sChrID  = list_sColumn[0]
            cVCF.nChrID  = process_chromo_ID(list_sColumn[0])
        else:
            if list_sColumn[0] == 'MT': continue
            if list_sColumn[0].startswith('<GL00'): continue
            if list_sColumn[0].startswith('.'): continue

            dict_sChrKey = {'X':'23', 'Y':'24'}
            cVCF.sChrID  = 'chr%s' % list_sColumn[0]

            if list_sColumn[0] in ['X', 'Y']:
                cVCF.nChrID  = int(dict_sChrKey[list_sColumn[0]])
            else:
                cVCF.nChrID  = int(list_sColumn[0])
        #if END:
        try: cVCF.nPos           = int(list_sColumn[1])
        except ValueError: continue
        cVCF.sDBSNP_ID      = list_sColumn[2]
        cVCF.sRefNuc        = list_sColumn[3]
        cVCF.sAltNuc        = list_sColumn[4]
        cVCF.fQual          = float(list_sColumn[5]) if list_sColumn[5] != '.' else list_sColumn[5]
        cVCF.sFilter        = list_sColumn[6]
        cVCF.sInfo          = list_sColumn[7]
        cVCF.sFormat        = list_sColumn[8]
        cVCF.list_sSamples  = list_sColumn[9:]

        dict_sVCFInfo       = dict(sInfo.split('=') for sInfo in cVCF.sInfo.split(';')[:-1] if len(sInfo.split('=')) == 2)
        cVCF.sLocale        = dict_sGENIC_alt[dict_sVCFInfo['Func.refGene'].split('\\x3b')[0]]
        cVCF.sSNVType       = dict_sVCFInfo['ExonicFunc.refGene'].split('\\x3b')[0]
        cVCF.sGeneSym       = dict_sVCFInfo['Gene.refGene'].split('\\x3b')[0]

        if sSEQDATA == 'WXS':
            if cVCF.sFilter != 'PASS': continue
        else:
            if 'PASS' not in sReadLine: continue
            if list_sColumn[2] == 'PASS':
                continue
        if bClustered:
            cVCF.nClusterID  = list_sColumn[-1]

        list_sOutput.append(cVCF)
    #loop END: sReadLine
    InFile.close()

    return list_sOutput
#def END: cVCF_parse_vcf_files
## endregion

##region class cClusterData

class cClusterData:
    def __init__(self):
        self.sChrID         = ''
        self.nPos           = 0
        self.sDBSNP_ID      = ''
        self.sRefNuc        = ''
        self.sAltNuc        = ''
        self.nVarCnt        = 0
        self.sInfo          = ''
        self.sPatID         = ''
        self.nClusterID     = 0
        self.sLocale        = ''
        self.sGeneSym       = ''
        self.fBackground    = 0.0
    #def END: __init__

def cClusterData_load_clustered_VCFs (sInFile):
    print_start_msg('Loading Clustered VCF', sInFile)
    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:
        # File Format - Altered VCF File Format
        # Column Number:     | 0       | 1         | 2          | 3        | 4
        # Column Description:| sChrID  | nPos      | sDBSNP_ID  | sRefNuc  | sAltNuc
        # Column Example:    | chr13   | 32906558  | rs79483201 | T        | A
        # Column Number:     | 5       | 6         | 7                     | 8
        # Column Description:| nVarCnt | sInfo     | sPatID                | nClusterID
        # Column Example:    | 4196    | GenicInfo | TCGA_LIHC_Mutect2_174 | 1

        if sReadLine.startswith('#'): continue  # SKIP Information Headers
        list_sColumn            = sReadLine.strip('\n').split('\t')

        cClust                = cClusterData()
        cClust.sChrID         = list_sColumn[0]
        cClust.nPos           = int(list_sColumn[1])
        cClust.sDBSNP_ID      = list_sColumn[2]
        cClust.sRefNuc        = list_sColumn[3]
        cClust.sAltNuc        = list_sColumn[4]
        cClust.nVarCnt        = int(list_sColumn[5])
        cClust.sInfo          = list_sColumn[6]
        cClust.sPatID         = list_sColumn[7]
        cClust.nClusterID     = int(list_sColumn[8])

        # Split the Info column by ';' then key:value by '='
        dict_sInfo      = dict(sInfo.split('=') for sInfo in cClust.sInfo.split(';')[:-1] if len(sInfo.split('=')) == 2)
        cClust.sLocale  = dict_sInfo['Func.refGene'].split('\\x3b')[0]
        cClust.sGeneSym = dict_sInfo['Gene.refGene'].split('\\x3b')[0].upper()
        cClust.sGeneID  = dict_sInfo['Gene.ensGene'].split('\\x3b')[0].upper()

        list_sOutput.append(cClust)
    #loop END: sReadLine

    #VS-Check
    if not list_sOutput:
        sys.exit('Invalid List : cClusterData_load_clustered_VCFs : list_sOutput size= %d' % len(list_sOutput))

    print_done_msg('Loading Clustered VCF', sInFile)
    return list_sOutput
#def END: cClusterData_load_clustered_VCFs



##endregion

## region class cRBP_Cluster

class cRBP_Cluster:pass

## endregion

##region class cHotspotData

class cHotspotData:

    def __init__(self):
        self.sChrID         = ''
        self.sStrand        = ''
        self.nStartPos      = 0
        self.nEndPos        = 0
        self.nSize          = 0
        self.nVarCnt        = 0
        self.nTotVarCnt     = 0
        self.fPvalue        = 0.0
        self.fQvalue        = 0.0
        self.sLocale        = ''
        self.fBackProb      = 0.0 # Sample background probability nTotVarCnt / genome (3 billion bps)

        self.list_sPatIDs   = []
        self.list_sCovered  = []
        self.list_cClusters = []
    #def END: __init__


def cHotspot_determine_back_probability(sVCFDir, list_cCluster, sChrID, nBackWindow, nStartPos, nEndPos):

    list_sOutput = []
    list_sPatID  = list(set([cClust.sPatID for cClust in list_cCluster]))

    for sPatID in list_sPatID:
        # Determine variant count in determined window
        sVCFFile           = '%s/%s.vcf.gz' % (sVCFDir, sPatID)
        nFlankSize         = int(nBackWindow / 2)

        if nStartPos > nFlankSize:
            nWinStart = nStartPos - nFlankSize - 1
        else: nWinStart = 1
        nWinEnd            = nEndPos + nFlankSize + 1

        ## Left Flank
        sScript            = 'tabix %s %s:%s-%s' % (sVCFFile, sChrID, nWinStart, nStartPos)
        sStdOut            = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout
        nVarCnt_left       = len([sReadLine for sReadLine in sStdOut])

        ## Right Flank
        sScript            = 'tabix %s %s:%s-%s' % (sVCFFile, sChrID, nEndPos, nWinEnd)
        sStdOut            = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout
        nVarCnt_right      = len([sReadLine for sReadLine in sStdOut])
        fBackground        = (nVarCnt_left + nVarCnt_right) / nBackWindow
        list_sOutput.append(fBackground)
    #loop END: sPatID

    #VS-Check
    if not list_sOutput:
        sys.exit('Invalid list : cHotspot_determine_back_probability : list_sOutput size= %d' % (len(list_sOutput)))

    return np.mean(list_sOutput)
# def END: cHotspot_determine_back_probability

## endregion

## region Class cCosmic
class cCOSMIC:
    def __init__(self):
        self.sGeneName   = ''
        self.sAccID      = ''
        self.nCDSLen     = 0
        self.sHGCNID     = ''   # SKIP for now
        self.sSample     = ''   # SKIP for now
        self.sSampleID   = ''   # SKIP for now
        self.sTumorID    = ''   # SKIP for now
        self.sPriSite    = ''   # primary site  ex) pancreas
        self.sSiteSub1   = ''   # SKIP for now
        self.sSiteSub2   = ''   # SKIP for now
        self.sSiteSub3   = ''   # SKIP for now
        self.sPriHist    = ''   # primary histology
        self.sHistSub1   = ''   # SKIP for now
        self.sHistSub2   = ''   # SKIP for now
        self.sHistSub3   = ''   # SKIP for now
        self.bGenomeWide = ''   # ex) y or n
        self.sMutaID     = ''   # SKIP for now
        self.sAltType    = ''   # ex) c.35G>T
        self.sAAType     = ''   # ex) p.G12V
        self.sMutaDescri = ''   # ex) Substitution - Missense
        self.sMutaZygo   = ''   # SKIP for now
        self.bLOH        = ''   # loss of heterzygosity ex) y or n
        self.sGRCh       = ''
        self.sChrom      = ''
        self.sPos        = 0
        self.sStrand     = ''
        self.bSNP        = ''   # ex) y and n
        self.sDelete     = ''   # ex) PATHOGENIC
    #def END : __init__

def cCos_parse_cosmic_consensus (sInFile):

    print_start_msg('Parsing COSMIC File', sInFile)

    dict_Check   = {}
    list_sOutput = []
    InFile       = open(sInFile, 'r')

    for i, sReadLine in enumerate(InFile):

        if sReadLine.startswith('Gene'):continue

        list_sColumn = sReadLine.strip('\n').split('\t')
        '''
        if i == 0:
            list_sHeader = list_sColumn
        elif i == 1:
            for i,(a,b) in enumerate(zip(list_sHeader, list_sColumn)):
                print('%s\t%s\t%s' % (i,a,b))
        else: break

        '''
        cCos             = cCOSMIC()
        cCos.sGeneName   = list_sColumn[0].upper()
        cCos.sAccID      = list_sColumn[1]
        cCos.nCDSLen     = int(list_sColumn[2])
        cCos.sHGCNID     = list_sColumn[3]
        cCos.sSample     = list_sColumn[4]
        cCos.sSampleID   = list_sColumn[5]
        cCos.sTumorID    = list_sColumn[6]
        cCos.sPriSite    = list_sColumn[7]
        cCos.sSiteSub1   = list_sColumn[8]
        cCos.sSiteSub2   = list_sColumn[9]
        cCos.sSiteSub3   = list_sColumn[10]
        cCos.sPriHist    = list_sColumn[11]
        cCos.sHistSub1   = list_sColumn[12]
        cCos.sHistSub2   = list_sColumn[13]
        cCos.sHistSub3   = list_sColumn[14]
        cCos.bGenomeWide = True if list_sColumn[15] == 'y' else False
        cCos.sMutaID     = list_sColumn[16]
        cCos.sAltType    = list_sColumn[17]
        cCos.sAAType     = list_sColumn[18]
        cCos.sMutaDescri = list_sColumn[19]
        cCos.sMutaZygo   = list_sColumn[20]
        cCos.bLOH        = True if list_sColumn[21] == 'y' else False
        cCos.sGRCh       = list_sColumn[22]
        cCos.sChrom      = 'chr%s' % list_sColumn[23].split(':')[0]
        list_sPosCheck   = list(set(list_sColumn[23].split(':')[1].split('-')))
        cCos.sDelete     = list_sColumn[26]
        cCos.sMutaStatus = list_sColumn[28]


        if len(list_sPosCheck) > 1:
            cCos.sPos    = list_sPosCheck[0]
        else:
            cCos.sPos    = ''.join(list_sPosCheck)

        cCos.sStrand     = list_sColumn[24]
        cCos.bSNP        = True if list_sColumn[25] == 'y' else False

        if cCos.sPriSite not in dict_Check:
            dict_Check[cCos.sPriSite] = 0
        dict_Check[cCos.sPriSite] += 1

        cCos.sDelete     = list_sColumn[26]

        if cCos.sDelete == '':
            cCos.sDelete = 'NotDefined'

        list_sOutput.append(cCos)

    #loop END: i, sReadLine

    InFile.close()

    #V-S Check:
    if not list_sOutput:
        sys.exit('Empty List : cCos_parse_cosmic_consensus : list_sOutput : Size = %d' % (len(list_sOutput)))

    return list_sOutput
#def END: cCos_parse_cosmic_consensus
#class END: cCosmic
## endregion

#region Class cMirBaseData
class cMirBaseData:
    def __init__(self):
        self.sGenome    = 'NULL'
        self.sMirName   = 'NULL'
        self.sChrID     = 'NULL'
        self.sStrand    = 'NULL'
        self.sFeature   = 'NULL'
        self.nStartPos  = 0
        self.nEndPos    = 0
    #def END: __init__
#class END: cMirBaseData

def cMirBaseData_read_GFF3_file (sInFile):

    cGenome    = cFasta(dict_sGENOME[sGENOME]) # For fetching primary miRNA sequence from genome
    list_cMir  = []
    InFile     = open(sInFile, 'r')
    for sReadLine in InFile:
        '''
        # GFF3 File Format
        # Column Number:        0             | 1       | 2                         | 3            | 4            | 5             |
        # Column Description:   contig        | SPACE   | feature                   | start pos    | end pos      | SPACE         |
        # Column Example:       chr1          | .       | miRNA_primary_transcript  | 30366        | 30503        | .             |
        #
        # Column Number:        6             | 7       | 8                                                                       |
        # Column Description:   strand        | SPACE   | attributes                                                              |
        # Column Example:       +             | .       | ID=MI0006363_1;accession_number=MI0006363;Name=hsa-mir-1302-2           |
        '''
        cMir = cMirBaseData()

        if sReadLine.startswith('#'):
            if sReadLine.startswith('# genome-build-id'):
                cMir.sGenome = sReadLine.strip('\n').split()[2]
        else:
            list_sColumn = sReadLine.strip('\n').split('\t')
            # V-S Check: Column List Size
            if len(list_sColumn) != nGFF3_COLUMN_SIZE:
                sys.exit('ERROR: GFF3 Column List Size= %d' % len(list_sColumn))

            cMir.sChrID     = list_sColumn[0]
            cMir.sFeature   = list_sColumn[2]
            cMir.nStartPos  = int(list_sColumn[3])  # 1-based closed
            cMir.nEndPos    = int(list_sColumn[4])  # 1-based closed
            cMir.sStrand    = list_sColumn[6]

            sAttribute      = list_sColumn[8].split(';')
            cMir.sMirName   = sAttribute[2].split('=')[1].lower()
            cMir.sMirSeq    = cGenome.fetch(cMir.sChrID, cMir.nStartPos-1, cMir.nEndPos, cMir.sStrand).upper()

            if cMir.sFeature == 'miRNA_primary_transcript':continue


            list_cMir.append(cMir)
        #if END: sReadLine.startswith
    #loop END: sReadLine
    InFile.close()

    # V-S Check: Empty List Check
    if not list_cMir:
        sys.exit('ERROR: cMir List Size= %d' % len(list_cMir))

    return list_cMir
#class END: cMirBaseData_read_GFF3_file
## endregion

## region Class cMiRNA
class cMiRNA:
    def __init__(self):
        self.sMirName           = ''
        self.sMirSeq            = ''
        self.dict_sST           = {}
        self.list_c3UTRs_CST    = []
        self.list_c3UTRs_NST    = []
        self.list_c3UTRs_CDNST  = []
        self.list_c3UTRs_NoSite = []
## endregion

## endregion


def main():
    #Qsub conditions
    bTestRun       = True
    sQueue         = 'optiplex.q'
    sTempScript    = copy_temp_core_script(sBASE_DIR)
    nJobs          = 1000

    #Analysis Conditions
    fRecurrency    = 0.01
    nClustDist     = 10
    nBackWinSize   = 100000 # or genome
    nBufferSize    = 20
    nVarWindow     = 0
    bStrictPeaks   = 1
    bAllRBPs       = 0


    ## Analysis Settings Conditions
    print(cyan('**********Analysis Conditions**********',1))
    print(cyan('Study and SeqData\t%s-%s'             % (sSTUDY, sSEQDATA),1))
    print(cyan('CellLine\t\t%s'                       % sCELLS,1))
    print(cyan('PeakBufferSize\t\t%s'                 % nBufferSize,1))
    print(cyan('Sample Cnt\t\t%s'                     % nSAMPLE_CNT,1))
    print(cyan('Recurrency\t\t%s'                     % int(fRecurrency*100),1))
    print(cyan('TopExpressed\t\t%s%%'                 % int(fTOP_EXPR*100),1))
    print(cyan('VariantWin\t\t%s'                     % nVarWindow,1))
    print(cyan('StrictPeaks\t\t%s'                    % bStrictPeaks,1))
    print(cyan('bAllRBPs\t\t%s'                       % bAllRBPs,1))
    print(cyan('***************************************',1))

    ## region Phase 1
    preprocess       (nClustDist, nBufferSize, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun)

    ## Only Applicable to Whole Exome Seq.
    #capture_filter   (nClustDist, nBufferSize, sTempScript, sQueue, bTestRun)

    determine_enrich (nVarWindow, nClustDist, nBufferSize, fRecurrency, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun)
    ## endregion

    ## region Phase 2 RBPxRBP, Cancer by Cancer analysis
    #rbp_by_rbp_analysis (nVarWindow, nClustDist, nBufferSize, fRecurrency, bStrictPeaks, bAllRBPs)

    #for sStudy in list_sCANCERs:
    #top_RBP_examination (sStudy, nVarWindow, nClustDist, nBufferSize, fRecurrency, bStrictPeaks, bAllRBPs)


    #determine_target_genes (nClustDist, nBufferSize, fRecurrency)
    ## endregion


    ''' OLD
    #sFileTag       = 'clust%s_back%s_buff%s_update' % (nClustDist, nBackWinSize, nBufferSize)

    #PRE-STEP 1: Basic Statistics of RBP Data (eCLIP and iCLIP)#
    #basic_stats_rbp                (list_cClipData)
    #buffer_rbps                    (list_cClipData, nBufferSize)

    #PRE-STEP 2: Run ANNOVAR on the VCF data in order to obtain genic location and other info#
    #qsub_annovar_vcfs              (sQueue, bTestRun) # Annotate VCFs

    #PRE-STEP 2-1: Bastic Stats of VCF data with ANNOVAR annotation information
    #basic_stats_vcf                ()
    #basic_stats_vcf_v2             () # For Testing

    #STEP 1: Prepare VCFs for clustering
    #output_vcfdata_for_clustering    ()

    #STEP 2: Cluster using BEDtools
    #qsub_bedtools_cluster            (nClustDist, sQueue, bTestRun)
    #basic_stats_clusters             (nClustDist)

    #STEP 3: Statistically evaluate clusters for hotspots
    #qsub_evaluate_clusters           (nClustDist, nBackWinSize, sTempScript, sQueue, bTestRun, nJobs)

    #STEP 4: Perform multiple test correction (FDR)
    #qsub_multipletest                (nClustDist, nBackWinSize, sTempScript, sQueue, bTestRun, nJobs)

    #STEP 5: Survey statistically significant (P<0.05) hotspots
    #basic_stats_hotspots             (nClustDist, nBackWinSize)

    #STEP 6: Overlap with Bedtools
    #qsub_get_rbp_variants            (list_cClipData, nClustDist, nBackWinSize, sTempScript, sQueue, bTestRun)

    #STEP 7-1: Basic Statistics of Clusters-RBPs P<0.05
    #basic_stats_clusters_rbps        (list_cClipData, nClustDist, nBackWinSize, fRecurrency)

    #STEP 7: Compile RBP data
    #compile_rbp_data                 (list_cClipData, nClustDist, nBackWinSize, sFileTag)

    #STEP 8: Compare miRNA target sites and (non-overlapping but in X-bp window)
    #determine_miRNA_targetsites      (sGenome, nClustDist, nBackWinSize, nBufferSize, sFileTag)

    #STEP 9-1: Overlap with Bedtools Part 2
    #output_rbpdata_forbuffer         (nClustDist, nBackWinSize, nBufferSize, sFileTag)

    #STEP 9-2: Overlap with Bedtools Part 2
    #qsub_get_rbp_variants_pt2        (list_cClipData, nClustDist, nBackWinSize, nBufferSize, sFileTag, sTempScript, sQueue, bTestRun)

    #STEP 9-3: Filter RBP-Hotspots that overlap with capture region
    #qsub_determine_rbp_coverage      (list_cClipData, nClustDist, nBackWinSize, nBufferSize, sFileTag, sTempScript, sQueue, bTestRun)

    #STEP 10: Determine RBP coverage
    #determine_rbp_coverage           (list_cClipData, nClustDist, nBackWinSize, nBufferSize, sFileTag)

    #STEP 11: Determine Enrichment by Genic Region
    #normalize_enrichment             (nClustDist, nBackWinSize, nBufferSize, fRecurrency, sFileTag)

    #STEP 8: Expression Data Comparison
    #evalute_expression_data (nClustDist, nBackWinSize)
    #evalute_expression_data_v2 (nClustDist, nBackWinSize, fRecurrency)

    '''
#def END: main

## region RBP-TCGA Phase 1

## region Preprocess

def preprocess_old (nClustDist, nBufferSize, bStrictPeaks, sTempScript, sQueue, bTestRun):
    # Proofread and Confirmed 1/21/17

    print(blue('Pre-process Data', 1))
    sPrePDir = '%s/01_preprocess_data'    % sBASE_DIR
    sWorkDir = '%s/%s'                    % (sPrePDir, sSTUDY)
    sClipDir = '%s/clipdata'              % sPrePDir
    sVCFDir  = '%s/VCFs-%s'               % (sWorkDir, sSEQDATA)
    sMirDir  = '%s/mir'                   % sWorkDir
    os.makedirs(sWorkDir, exist_ok=True)
    os.makedirs(sClipDir, exist_ok=True)
    os.makedirs(sVCFDir, exist_ok=True)
    os.makedirs(sMirDir, exist_ok=True)

    '''Process RefSeq Data'''
    sRefFile     = '%s/refseq.dat'    % sPrePDir
    dict_cRefSeq = process_refseq_data (sPrePDir, sRefFile)

    '''Process Expr Data''' # For checking RBP expression
    #process_expr_data (sRefFile)

    '''Process eCLIP Data'''
    #qsub_rename_samples (sRefFile, sTempScript, sQueue, bTestRun)
    #process_eclip_data  (sWorkDir, sClipDir, sRefFile, bStrictPeaks, sTempScript, sQueue, bTestRun) # Confirmed 1/21/17

    '''Process VCF Data'''
    #process_vcf_data    (sVCFDir, sRefFile, nClustDist, sTempScript, sQueue, bTestRun)

    '''Preprocess TCGA FPKM files'''
    #qsub_preprocess_fpkm (sTempScript, sQueue, bTestRun)

    '''Determine miRNA target sites'''
    #process_mir_data    (sMirDir, sRefFile) # Confirmed 1/23/17
#def END: preprocess

def preprocess (nClustDist, nBufferSize, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun):
    # Proofread 1/23/17
    print(blue('Pre-process Data', 1))

    sWorkDir = '%s/01_preprocess_data' % sBASE_DIR
    os.makedirs(sWorkDir, exist_ok=True)

    ### Process RefSeq Data ###
    sRefFile     = '%s/refseq.dat'     % sWorkDir
    #dict_cRefSeq = process_refseq_data (sWorkDir, sRefFile)
    ###########################

    ### Process Expr Data #####
    #process_expr_data (sRefFile)
    ###########################

    ### Process eCLIP Data ####
    sClipDir    = '%s/clipdata'        % sWorkDir
    os.makedirs(sClipDir, exist_ok=True)
    #qsub_rename_samples (sRefFile, sTempScript, sQueue, bTestRun) # For renaming gene names bedfiles
    #process_eclip_data  (sClipDir, sRefFile, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun)
    ###########################

    ### Process VCF Data #####
    sVCFDir    = '%s/%s/VCFs-%s'           % (sWorkDir, sSTUDY, sSEQDATA)
    os.makedirs(sVCFDir, exist_ok=True)
    process_vcf_data    (sVCFDir, nClustDist, sTempScript, sQueue, bTestRun)
    ##########################

    ### Preprocess TCGA FPKM files ###
    #qsub_preprocess_fpkm (sTempScript, sQueue, bTestRun)
    ##################################

    ### Determine miRNA target sites ###
    sMirDir    = '%s/mir'               % sWorkDir
    os.makedirs(sMirDir, exist_ok=True)
    #process_mir_data    (sMirDir, sRefFile)
    ####################################
#def END: preprocess



## region RefSeq

def process_refseq_data (sOutDir, sOutFile):
    print(green('Process RefSeq data', 1))

    dict_cRefSeq = cRef_parse_refflat_line(sREFSEQFILE, sOutFile)

    print('Total RefSeq Genes', len(dict_cRefSeq))
    OutFile = open('%s/refseqcheck.txt' % sOutDir, 'w')
    for sGeneSym in dict_cRefSeq:
        cRef = dict_cRefSeq[sGeneSym]
        #Out format: GeneSym NMID 5UTRsize ORFsize 3UTRsize
        sOut = '%s\t%s\t%s\t%s\t%s\n' % (sGeneSym, cRef.sNMID, len(cRef.s5UTRSeq),
                                         len(cRef.sORFSeq), len(cRef.s3UTRSeq))
        OutFile.write(sOut)
    #loop END: sGeneSym
    OutFile.close()

    return dict_cRefSeq
#def process_refseq_data

## endregion

## region Expr

def process_expr_data (sRefFile):

    print(green('Process Expr data', 1))

    ## Using transcriptome data in K562 cells ##
    #reformat_encode_gtf_file (sRefFile)
    ############################################

    dict_sExpr     = {'K562':  load_expr_data('K562'),
                      'HepG2': load_expr_data ('HepG2')}
    list_cClipData = get_clipdata_list(1)

    dict_sRBPID    = {}
    for cClipData in list_cClipData:

        sKey = cClipData.sGeneSym.split('-')[0]

        if sKey not in dict_sRBPID:
            dict_sRBPID[sKey] = []
        dict_sRBPID[sKey].append(cClipData)
    #loop END: cClipData

    sOutFile = '%s/RBP_Expression.txt' % sBASE_DIR
    OutFile  = open(sOutFile, 'w')

    sHeader  = '%s\t%s\t%s\t%s\n' % ('RBP', 'Cells', 'K562', 'HepG2')
    OutFile.write(sHeader)

    for sRBPID in dict_sRBPID:

        if len(dict_sRBPID[sRBPID]) == 4:

            nTotal1 = len(dict_sExpr['K562'])
            nTotal2 = len(dict_sExpr['HepG2'])
            fRank1, fRPKM1 = dict_sExpr['K562'][sRBPID]
            fRank2, fRPKM2 = dict_sExpr['HepG2'][sRBPID]
            sOut    = '%s\t%s\t%s\t%s\t%s\t%s\n' \
                      % (sRBPID, 'Both', fRank1/nTotal1, fRank2/nTotal2, fRPKM1, fRPKM2)
            OutFile.write(sOut)

        else:
            sCells = [cClipData.sCellLine for cClipData in dict_sRBPID[sRBPID]][0]

            nTotal = len(dict_sExpr[sCells])
            fRank, fRPKM = dict_sExpr[sCells][sRBPID]
            sOut = '%s\t%s\t%s\t%s\n' % (sRBPID, sCells, fRank/nTotal, fRPKM)
            OutFile.write(sOut)
        #if END:
    #loop END: sRBPID
    OutFile.close()
#def END: process_expr_data


def reformat_encode_gtf_file (sRefFile):
    #Reformat transcriptome data for consistency
    sWorkDir = '%s/exprdata' % sBASE_DIR
    for i in [1, 2]:
        sInFile      = '%s/K562_rep%s.gene.gtf' % (sWorkDir, i)
        dict_sOutput = load_gtf_file(sInFile)

        sOutFile     = '%s/K562_rep%s.rpkm.gtf' % (sWorkDir, i)
        list_sOutput = [[dict_sOutput[sGeneID][0], sGeneID, dict_sOutput[sGeneID][1]]
                        for sGeneID in dict_sOutput]
        list_sOutput = sorted(list_sOutput, key=lambda  e:e[2], reverse=True)

        OutFile      = open(sOutFile, 'w')
        for sGeneSym, sGeneID, fFPKM in list_sOutput:

            sOut = '%s\t%s\t%s\t%s\n' % (sGeneSym, sGeneID, 'filler', fFPKM)
            OutFile.write(sOut)
        #loop END: sGeneSym, sGeneID, fFPKM
        OutFile.close()
    # loop END: i
#def END: reformat_encode_gtf_files


def load_gtf_file (sInFile):

    dict_sOutput = {}
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:

        if sReadLine.startswith('##'): continue #Skip headers

        ''' File Format
        0 sChrID          chr1
        1 sAnnotation     HAVANA
        2 sFeature        gene
        3 nStartPos       11869   *1-based
        4 nEndPos         14409
        5 sScore          .       *Not used
        6 sStrand         +
        7 nGenomicPhase   .
        8 sAddInfo        gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1";  level 2; havana_gene "OTTHUMG00000000961.2";
        '''
        list_sColumn    = sReadLine.strip('\n').split('\t')

        sChrID          = list_sColumn[0]
        sAnnotation     = list_sColumn[1]
        sFeature        = list_sColumn[2]

        if sFeature != 'gene': continue # Skip other redundant info (exon, CDS, transcript, UTR...)

        nStartPos       = int(list_sColumn[3])
        nEndPos         = int(list_sColumn[4])
        sScore          = list_sColumn[5]
        sStrand         = list_sColumn[6]
        nGenomicPhase   = list_sColumn[7]
        sGeneInfo       = list_sColumn[8].replace('"', '')

        dict_sGeneInfo  = dict([list(filter(None, sInfo.split(' '))) for sInfo in sGeneInfo.split(';') if list(filter(None, sInfo.split(' ')))])

        sGenesym        = dict_sGeneInfo['gene_name'].upper()
        sGeneID         = dict_sGeneInfo['gene_id'].upper()
        fFPKM           = float(dict_sGeneInfo['FPKM'])

        if sGeneID not in dict_sOutput:
            dict_sOutput[sGeneID] = ''
        dict_sOutput[sGeneID] = [sGenesym, fFPKM]
    #loop END: sReadLine

    #V-S Check:
    if not dict_sOutput:
        sys.exit('Invalid Dictionary : load_gtf_file : dict_sOutput size= %d' % len(dict_sOutput))

    return dict_sOutput
#def END: load_gtf_file


def load_expr_data (sCells, fTop50=0):

    sWorkDir     = '%s/exprdata' % sBASE_DIR
    dict_fExpr   = {}

    if sCells == 'K562': # Has replicates
        for i in [1, 2]:
            sInFile    = '%s/%s_rep%s.rpkm.gtf' % (sWorkDir, sCells, i)
            InFile     = open(sInFile, 'r')
            list_fExpr = [list(sReadLine.strip('\n').split('\t')) for sReadLine in InFile]
            InFile.close()

            for sGeneSym, sGeneID, sFiller, fFPKM in list_fExpr:
                if sGeneSym not in dict_fExpr:
                    dict_fExpr[sGeneSym] = []
                dict_fExpr[sGeneSym].append(float(fFPKM))
            #loop END: sGeneSym, sGeneID, sFiller, fFPKM
        #loop END: i

    elif sCells == 'HepG2':

        sInFile    = '%s/%s_rep1.rpkm.gtf' % (sWorkDir, sCells)
        InFile     = open(sInFile, 'r')
        list_fExpr = [list(sReadLine.strip('\n').split('\t')) for sReadLine in InFile]
        InFile.close()

        for sGeneSym, sGeneID, sFiller, fFPKM in list_fExpr:
            if sGeneSym not in dict_fExpr:
                dict_fExpr[sGeneSym] = []
            dict_fExpr[sGeneSym].append(float(fFPKM))
        #loop END: sGeneSym, sGeneID, sFiller, fFPKM
    #if END: sCells

    list_sOutput = [[sGeneSym, np.mean(dict_fExpr[sGeneSym])] for sGeneSym in dict_fExpr]
    list_sOutput = sorted(list_sOutput, key=lambda e:e[1], reverse=True)
    if fTop50:
        list_sOutput = list_sOutput[:int(len(list_sOutput) * 0.5)]

    #print(len(list_sOutput))

    dict_sOutput = {}
    for i, (sGeneSym, fFPKM) in enumerate(list_sOutput):

        if sGeneSym not in dict_sOutput:
            dict_sOutput[sGeneSym] = []
        dict_sOutput[sGeneSym] = [i, fFPKM]
    #loop END: i, (sGeneSym, fFPKM)

    return dict_sOutput
#def END: list_sOutput

## endregion



## region eCLIP

def qsub_rename_samples (sRefFile, sTempScript, sQueue, bTestRun):

    print(green('Process eCLIP data', 1))

    '''Load CLIP data'''
    sClipBedDir    = '%s/sorted_bed'                 % sECLIP_DIR
    sOutDir        = '%s/sorted_bed_renamed'         % sECLIP_DIR
    sJobName       = 'Jinman.FixNames.eCLIP.%s'      % sCELLS
    sLogDir        = '%s/log/%s/%s'                  % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    list_cClipData = get_clipdata_list()
    dict_cClipData = {}

    for cClipData in list_cClipData:

        sKey = cClipData.sGeneSym.split('-')[0]

        if sKey not in dict_cClipData:
            dict_cClipData[sKey] = []
        dict_cClipData[sKey].append(cClipData.sAccID)
    #loop END: cClipData

    for sRBPID in dict_cClipData:
        sScript        = '%s rename_samples %s %s %s %s ' % \
                         (sTempScript, sRBPID, sRefFile, sClipBedDir, sOutDir)
        if bTestRun: print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s'
                      % (sScript, sLogDir, sQueue, sJobName, sRBPID))
        #if END
    #loop END: cClipData
#def END: qsub_rename_samples


def rename_samples (sRBPID, sRefFile, sInDir, sOutDir):

    print('Opening Pickled', sRefFile)
    OutFile        = open(sRefFile, 'rb')
    dict_cRefSeq   = pickle.load(OutFile)
    OutFile.close()

    list_cClipData = get_clipdata_list()
    dict_cClipData = {}

    for cClipData in list_cClipData:

        sKey = cClipData.sGeneSym.split('-')[0]

        if sKey not in dict_cClipData:
            dict_cClipData[sKey] = {}

        if cClipData.sCellLine not in dict_cClipData[sKey]:
            dict_cClipData[sKey][cClipData.sCellLine] = []

        dict_cClipData[sKey][cClipData.sCellLine].append(cClipData)
    #loop END: cClipData

    for sCellLine in dict_cClipData[sRBPID]:

        for i, cClipData in enumerate(dict_cClipData[sRBPID][sCellLine]):

            sFileID         = cClipData.sAccID

            sSampleKey      = '%s_%s_rep0%s' % (sRBPID, sCellLine, i+1)

            sBedFile        = '%s/%s.bed'    % (sInDir, sFileID)
            sOutBed         = '%s/%s.bed'    % (sOutDir, sFileID)

            list_cPeakData  = cPeakData_parse_bedfile(sBedFile, dict_cRefSeq, 1)

            print(sOutBed, len(list_cPeakData))

            OutFile         = open(sOutBed, 'w')
            for cPeak in list_cPeakData:
                sOut = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                       (cPeak.sChrID, cPeak.nStartPos, cPeak.nEndPos, sSampleKey,
                        cPeak.nScore, cPeak.sStrand, cPeak.fSignal, cPeak.fPvalue,
                        cPeak.sNotUsed1, cPeak.sNotUsed2)
                OutFile.write(sOut)
            #loop END: cPeak
            OutFile.close()
        #loop END: i, cClipData
    #loop END: sCellLine
#def END: rename_samples


def process_eclip_data (sClipDir, sRefFile, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun):
    print(green('Process eCLIP data', 1))

    sClipBedDir      = '%s/sorted_bed'    % sECLIP_DIR
    list_cClipData   = get_clipdata_list(bAllRBPs)

    ## For quick analysis of molecular pathways involved with the RBPs
    #list_sTargetRBPs = pathway_survey (list_cClipData)
    #list_cClipData   = [cClipData for cClipData in list_cClipData if cClipData.sAccID in list_sTargetRBPs]
    #print(len(list_sTargetRBPs))
    #print(len(list_cClipData))

    '''Process eCLIP data'''
    qsub_process_eclipdata (sClipDir, sRefFile, list_cClipData, sClipBedDir, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun)

    '''Basic Stats of eCLIP data'''
    get_peak_distribution (sClipDir, list_cClipData, bStrictPeaks, 0)
#def END: process_eclip_data


def qsub_process_eclipdata (sClipDir, sRefFile, list_cClipData, sClipBedDir, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun):

    sPeakCall       = 'strict' if bStrictPeaks else 'allpeaks'
    sRBPGroup       = 'AllRBPs' if bAllRBPs else sCELLS
    sJobName        = 'Jinman.ProcessFilter.eCLIP.%s.%s' % (sRBPGroup, sPeakCall)
    sLogDir         = '%s/log/%s/%s'                     % (sBASE_DIR, sJobName, sTIME_STAMP)
    sOutDir         = '%s/%s/%s-peaks'                   % (sClipDir, sPeakCall, sRBPGroup)
    sOutbyChr       = '%s/%s/%s-peaks_bychr'             % (sClipDir, sPeakCall, sRBPGroup)

    os.makedirs(sLogDir, exist_ok=True)
    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sOutbyChr, exist_ok=True)

    list_sBedFiles  = []

    for cClipData in list_cClipData:
        sBedFile       = '%s/%s.bed'                           % (sClipBedDir, cClipData.sAccID)
        list_sBedFiles.append(sBedFile)

        sOutFile       = '%s/%s.peakdata'                      % (sOutDir, cClipData.sAccID)
        sScript        = '%s process_eclipdata %s %s %s %s %s %s'      \
                         % (sTempScript, cClipData.sAccID, sRefFile, sBedFile, sOutFile, sOutbyChr, bStrictPeaks)
        if bTestRun: print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s'
                      % (sScript, sLogDir, sQueue, sJobName, cClipData.sAccID))
        #if END
    #loop END: cClipData

    sJobName       = 'Jinman.%s.%s.Bedtools.Merge.%s' % (sSTUDY, sSEQDATA, sPeakCall)
    sLogDir        = '%s/log/%s/%s'                   % (sBASE_DIR, sJobName, sTIME_STAMP)
    sTempDir       = '%s/temp'                        % sClipDir
    os.makedirs(sTempDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    sOutFile       = '%s/%s/%s_combined.bed'               % (sClipDir, sPeakCall, sRBPGroup)
    sOutFile2      = '%s/%s/%s_combined.sorted.bed'        % (sClipDir, sPeakCall, sRBPGroup)
    sOutFile3      = '%s/%s/%s_combined.sorted.merged.bed' % (sClipDir, sPeakCall, sRBPGroup)

    if bStrictPeaks: sKey = 200 #Filter out, Keep score == 1000
    else:            sKey = 'None'

    sScript        = 'cat %s | awk \'{ if (\\$5 != \\\"%s\\\") print }\' > %s ;'\
                                                        % (' '.join(list_sBedFiles), sKey, sOutFile)
    sScript       += 'sort -k1,1 -k2,2n -T %s %s > %s;' % (sTempDir, sOutFile, sOutFile2)
    sScript       += 'bedtools2 merge -c 4,5,6,7,8,9,10 -o distinct  -i %s > %s;' \
                                                        % (sOutFile2, sOutFile3)
    sScript       += '%s parse_merged_bed %s %s;'       % (sTempScript, sOutFile3, sRefFile)

    if bTestRun: print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s'
                      % (sScript, sLogDir, '4_730.q', sJobName))
    #if END:
#def END: qsub_process_eclipdata


def process_eclipdata (sAccID, sRefFile, sInFile, sOutFile, sOutbyChr, bStrictPeaks):

    print('Opening Pickled', sRefFile)

    bStrictPeaks    = int(bStrictPeaks)
    InFile          = open(sRefFile, 'rb')
    dict_cRefSeq    = pickle.load(InFile)
    InFile.close()

    dict_sChrID     = {sChrID:[] for sChrID in list_sCHRIDs}
    list_cPeakData  = cPeakData_parse_bedfile(sInFile, dict_cRefSeq)

    if bStrictPeaks: # Peak Filter
        list_cPeakData = [cPeak for cPeak in list_cPeakData if cPeak.nScore == 1000]

    print(sInFile, len(list_cPeakData))

    OutFile         = open(sOutFile, 'wb')
    pickle.dump(list_cPeakData, OutFile)
    OutFile.close()

    for cPeak in list_cPeakData:

        if bStrictPeaks: # Peak Filter
            if cPeak.nScore != 1000: continue

        if cPeak.sChrID not in dict_sChrID:
            dict_sChrID[cPeak.sChrID] = []
        dict_sChrID[cPeak.sChrID].append(cPeak)
    #loop END:

    for sChrID in list_sCHRIDs:
        print(sChrID, len(dict_sChrID[sChrID]))

        sOutChr = '%s/%s' % (sOutbyChr, sChrID)
        os.makedirs(sOutChr, exist_ok=True)

        OutFile         = open('%s/%s.%s.peakdata' % (sOutChr, sChrID, sAccID), 'wb')
        pickle.dump(dict_sChrID[sChrID], OutFile)
        OutFile.close()
    #loop END: sChrID

#def END: process_eclip_data


def parse_merged_bed (sInFile, sRefFile):

    print('Opening Pickled', sRefFile)
    OutFile         = open(sRefFile, 'rb')
    dict_cRefSeq    = pickle.load(OutFile)
    OutFile.close()

    sOutFile        = sInFile.replace('.bed', '.dat')
    list_cPeakData  = cPeakData_parse_bedfile_alt(sInFile, dict_cRefSeq)

    print(sInFile, len(list_cPeakData))
    OutFile         = open(sOutFile, 'wb')
    pickle.dump(list_cPeakData, OutFile)
    OutFile.close()
#def END: parse_merged_bed


def get_peak_distribution (sClipDir, list_cClipData, bStrictPeaks, bFiltered=0):
    print_start_msg('Basic Stats for RBP : Peak Locale')

    sPeakCall     = 'strict' if bStrictPeaks else 'allpeaks'
    dict_sExpr    = {'K562':  load_expr_data('K562', 1),
                     'HepG2': load_expr_data ('HepG2', 1)}
    dict_sLocale  = {}
    dict_nPeakCnt = {}
    sOutFile      = '%s/PeaksPerRBP%s%s.txt' % (sClipDir, '.cap' if bFiltered else '',
                                              '.strict' if bStrictPeaks else '.allpeaks')
    sPeakDir      = '%s/%s/peaks'            % (sClipDir, sPeakCall)

    OutFile       = open(sOutFile, 'w')
    list_sLocale  = ['exonic', '3UTR', '5UTR', 'intergenic', 'intronic']
    for cClipData in list_cClipData:

        # Load expression data
        dict_fFPKM     = dict_sExpr[cClipData.sCellLine]

        InFile         = open('%s/%s.peakdata' % (sPeakDir, cClipData.sAccID), 'rb')
        list_cPeakData = pickle.load(InFile)
        InFile.close ()

        OutFile.write('%s\t%s\n' % (cClipData.sGeneSym, len(list_cPeakData)))

        if cClipData.sGeneSym not in dict_nPeakCnt:
            dict_nPeakCnt[cClipData.sGeneSym] = {sLocale:0 for sLocale in list_sLocale}

        for cPeak in list_cPeakData:
            if cPeak.sLocale == '':
                cPeak.sLocale = 'intergenic'
            else:
                try: dict_fFPKM[cPeak.sGeneSym]
                except KeyError: continue
                pass

            dict_nPeakCnt[cClipData.sGeneSym][cPeak.sLocale] += 1

            if cPeak.sLocale not in dict_sLocale:
                dict_sLocale[cPeak.sLocale] = []
            dict_sLocale[cPeak.sLocale].append(cPeak)
        #loop END: cPeak
    #loop END: cClipData
    OutFile.close()

    sOutFile     = '%s/PeakDistPerRbp%s%s.txt' % (sClipDir, '.cap' if bFiltered else '',
                                              '.strict' if bStrictPeaks else '.allpeaks')
    OutFile2     = open(sOutFile, 'w')
    for sGeneSym in dict_nPeakCnt:
        list_nVarCnts = [dict_nPeakCnt[sGeneSym][sLocale] for sLocale in list_sLocale]
        sOut = '%s\t%s\n' % (sGeneSym, '\t'.join([str(nCnt) for nCnt in list_nVarCnts]))
        OutFile2.write(sOut)
    OutFile2.close()

    dict_sTargetGenes = {}
    for sLocale in list_sLocale:
        if 'UTR' not in sLocale: continue

        for cPeak in dict_sLocale[sLocale]:
            if cPeak.sGeneSym not in dict_sTargetGenes:
                dict_sTargetGenes[cPeak.sGeneSym] = 0
            dict_sTargetGenes[cPeak.sGeneSym] += 1
        #loop END: cPeak
    #loop END: sLocale
    print('UTRs', len(dict_sTargetGenes))

    #for sLocale in list_sLocale:
    #    print(sLocale, len(dict_sLocale[sLocale]))

    print_done_msg('Basic Stats for RBP : Peak Distribution')
#def END: get_peak_distribution

## endregion

## region VCF

def process_vcf_data (sWorkDir, nClustDist, sTempScript, sQueue, bTestRun):
    print(green('Process VCF data', 1))

    '''Run ANNOVAR on the VCF data in order to obtain genic location and other info'''
    #qsub_annovar_vcfs              (sWorkDir, sTempScript, sQueue, bTestRun) # Annotate VCFs

    ''''Basic Stats of VCF data with ANNOVAR annotation information'''
    #qsub_basic_stats_vcf           (sWorkDir, sTempScript, sQueue, bTestRun)  # Annotate VCFs
    #basic_stats_vcf                (sWorkDir)

    '''OLD Prepare VCFs for clustering'''
    #output_vcfdata_for_clustering   (sWorkDir)

    '''OLD Cluster using BEDtools'''
    #qsub_bedtools_cluster            (sWorkDir, nClustDist, sQueue, bTestRun)
#def END: process_vcf_data


def qsub_annovar_vcfs (sWorkDir, sTempScript, sQueue, bTestRun):
    print(green('Run ANNOVAR on the VCF data in order to obtain genic location and other info', 1))
    sVCFDir     = '%s/vardata/%s/%s/VCFs_hg19' % (sBASE_DIR, sSTUDY, sSEQDATA)
    sOutDir     = '%s/VCFs'                    % sWorkDir
    sChrIDDir   = '%s/VCFs_bychr'              % sWorkDir
    sTempDir    = '%s/VCFs/temp'               % sWorkDir
    sJobName    = 'Jinman.ANNOVAR.VCF.%s.%s'   % (sSEQDATA, sSTUDY)
    sLogDir     = '%s/log/%s/%s'               % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sChrIDDir, exist_ok=True)
    os.makedirs(sTempDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    # Load Patient IDs
    list_sPatIDs = get_patientIDs(sSTUDY, sSEQDATA)

    for sPatID in list_sPatIDs:
        sVCFFile  = '%s/%s.vcf.gz'             % (sVCFDir, sPatID)
        sAnnoFile = '%s/%s.hg19_multianno.vcf' % (sOutDir, sPatID)

        # Annotate with ANNOVAR
        if not os.path.isfile(sAnnoFile):
            sScript     = '%s %s '           % (sANNOVAR, sVCFFile)
            sScript    += '/extdata6/Jinman/util/annovar/humandb '
            sScript    += '-buildver hg19 '
            sScript    += '--vcfinput '
            sScript    += '--outfile %s/%s ' % (sOutDir, sPatID)
            sScript    += '--tempdir %s '    % sTempDir
            sScript    += '-protocol refGene,ensGene '
            sScript    += '-operation  g,g, ;'
            sScript    += '%s sep_by_chr %s %s %s; ' % (sTempScript, sPatID, sAnnoFile, sChrIDDir)
        else:
            sScript     = '%s sep_by_chr %s %s %s; ' % (sTempScript, sPatID, sAnnoFile, sChrIDDir)
        #if END:

        sInFile  = '%s/%s.hg19_multianno.vcf'    % (sOutDir, sPatID)
        sZipFile = '%s/%s.hg19_multianno.vcf.gz' % (sOutDir, sPatID)

        if os.path.isfile(sZipFile): continue
        sScript += 'cat %s | egrep PASS | awk \'{ if (\\$2 != \\\"PASS\\\") print }\' |  bgzip -c > %s;' % (
        sInFile, sZipFile)
        sScript += 'tabix -p vcf %s;' % sZipFile

        if bTestRun:
            print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s'
                      % (sScript, sLogDir, sQueue, sJobName, sPatID))
        #if END:
    #loop END: sPatID
#def END: qsub_annovar_vcfs


def qsub_basic_stats_vcf (sWorkDir, sTempScript, sQueue, bTestRun):
    print(green('Run ANNOVAR on the VCF data in order to obtain genic location and other info', 1))
    sInDir     = '%s/VCFs'                     % sWorkDir
    sStatDir   = '%s/VCFs_basics'              % sWorkDir
    sPrevJob   = 'Jinman.ANNOVAR.VCF.%s'       % sSEQDATA
    sJobName   = 'Jinman.VCF.Stats.%s.%s'      % (sSTUDY, sSEQDATA)
    sLogDir    = '%s/log/%s/%s'                % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sStatDir, exist_ok=True)

    # Load Patient IDs
    list_sPatIDs  = get_patientIDs(sSTUDY, sSEQDATA)

    for sPatID in list_sPatIDs:

        sAnnoFile = '%s/%s.hg19_multianno.vcf' % (sInDir, sPatID)

        # Annotate with ANNOVAR
        sScript   = '%s basic_stats_vcf_indy %s %s %s' % (sTempScript, sPatID, sAnnoFile, sStatDir)

        if bTestRun:
            print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s -hold_jid %s.%s'
                      % (sScript, sLogDir, sQueue, sJobName, sPatID, sPrevJob, sPatID))
        #if END:
    #loop END: sPatID
#def END: qsub_annovar_vcfs


def sep_by_chr (sPatID, sInFile, sOutDir):
    print(green('Basic Stats of VCF data with ANNOVAR annotation information', 1))

    list_cVCF    = cVCF_parse_vcf_files (sInFile)
    dict_sChrID  = {}
    for cVCF in list_cVCF:
        sKey          = cVCF.sChrID
        if sKey not in dict_sChrID:
            dict_sChrID[sKey] = []
        dict_sChrID[sKey].append(cVCF)
    #loop END: sPatID

    for sChrID in dict_sChrID:

        sChrDir = '%s/%s' % (sOutDir, sChrID)
        os.makedirs(sChrDir, exist_ok=True)

        print(sChrID, len(dict_sChrID[sChrID]))
        sOutFile = '%s/%s.%s.hg19_multianno.dat' % (sChrDir, sPatID, sChrID)
        OutFile  = open(sOutFile, 'wb')
        pickle.dump(dict_sChrID[sChrID], OutFile)
        OutFile.close()
    #loop END: sChrID
#def END: basic_stats_vcf


def basic_stats_vcf_indy (sPatID, sInFile, sOutDir, bFiltered=0):
    print(green('Basic Stats of VCF data with ANNOVAR annotation information', 1))

    list_sLocale = ['exonic', '5UTR', '3UTR', 'intergenic', 'intronic']

    if bFiltered:
        InFile    = open(sInFile, 'rb')
        list_cVCF = pickle.load(InFile)
        InFile.close()
    else:
        list_cVCF = cVCF_parse_vcf_files (sInFile)
    #if END:
    print(sPatID, len(list_cVCF))

    dict_sVarDist = {sLocale:0 for sLocale in list_sLocale}
    list_sSNVType = [len(list_cVCF), len([cVCF for cVCF in list_cVCF
                                          if cVCF.sSNVType == 'nonsynonymous_SNV'])]
    for cVCF in list_cVCF:
        sKey = dict_sGENIC_alt[cVCF.sLocale]
        if sKey not in dict_sVarDist:
            dict_sVarDist[sKey] = 0
        dict_sVarDist[sKey] += 1
    #loop END: sPatID

    OutFile = open('%s/%s_stats%s.txt' % (sOutDir, sPatID, '_cap' if bFiltered else ''), 'w')

    #OutFile.write('#%s\ttotvar\tnonSynCnt\n' % sPatID)
    OutFile.write('nonSynCnt-%s\t%s\t%s\n' % (sPatID, list_sSNVType[0], list_sSNVType[1]))

    # VarDist by Patient
    #OutFile.write('#%s\t%s\n' % (sPatID, '\t'.join(list_sLocale)))
    #sOut= '%s\t%s\n' % (sPatID, '\t'.join([str(dict_sVarDist[sLocale]) for sLocale in list_sLocale]))
    #OutFile.write(sOut)
    OutFile.close()
#def END: basic_stats_vcf


def load_cRefSeq (sInFile):
    dict_sOutput = {}
    InFile = open(sInFile, 'r')
    for sReadLine in InFile:
        list_sColumn  = sReadLine.strip('\n').split('\t')
        cRef          = cRefSeq()
        cRef.sGeneSym = list_sColumn[0].upper()
        cRef.sNMID    = list_sColumn[1]

        if cRef.sNMID not in dict_sOutput:
            dict_sOutput[cRef.sGeneSym] = ''
        dict_sOutput[cRef.sGeneSym] = cRef.sNMID
    #loop END: sReadLine

    return dict_sOutput
#def END: load_cRefSeq


def basic_stats_vcf (sWorkDir, bFiltered=0):
    print(green('Basic Stats of VCF data with ANNOVAR annotation information', 1))

    '''Deteremine regional sizes for normalization '''
    sVCFDir         = '%s/VCFs'      % sWorkDir
    sRefSeqFile     = '/extdata6/Jinman/reference_genome/091613_4WT.txt'
    dict_sRefSeq    = load_cRefSeq(sRefSeqFile)

    # Load Patient IDs
    list_sPatIDs  = get_patientIDs(sSTUDY, sSEQDATA)
    #
    dict_sLocale  = {}
    nTotalVar     = 0
    #sOutFile      = '%s/VariantsPerPatients%s.txt' % (sWorkDir, '_cap' if bFiltered else '')
    #OutFile       = open(sOutFile, 'w')
    dict_sVarDist = {}
    list_tLOD     = []
    list_sLocale  = ['exonic', '5UTR', '3UTR', 'intergenic', 'intronic']

    dict_sSNVType = {}
    dict_sGeneSym = {}
    for sPatID in list_sPatIDs:

        if sPatID not in dict_sSNVType:
            dict_sSNVType[sPatID] = []

        if bFiltered:
            sVCFFile  = '%s/%s.hg19_multianno.dat' % (sVCFDir, sPatID)
            InFile    = open(sVCFFile, 'rb')
            list_cVCF = pickle.load(InFile)
            list_cVCF = [cVCF for cVCF in list_cVCF if cVCF.sFilter == 'PASS']

            InFile.close()
        else:
            sVCFFile  = '%s/%s.hg19_multianno.vcf'  % (sVCFDir, sPatID)
            list_cVCF = cVCF_parse_vcf_files (sVCFFile)

            print(sPatID, len(list_cVCF))
        #if END:
    sys.exit()
    for sPatID in list_sPatIDs:

        if sPatID not in dict_sSNVType:
            dict_sSNVType[sPatID] = []

        if bFiltered:
            sVCFFile  = '%s/%s.hg19_multianno.dat' % (sVCFDir, sPatID)
            InFile    = open(sVCFFile, 'rb')
            list_cVCF = pickle.load(InFile)
            list_cVCF = [cVCF for cVCF in list_cVCF if cVCF.sFilter == 'PASS']

            InFile.close()
        else:
            sVCFFile  = '%s/%s.hg19_multianno.vcf'  % (sVCFDir, sPatID)
            list_cVCF = cVCF_parse_vcf_files (sVCFFile)

            print(sPatID, len(list_cVCF))
        #if END:

        nTotalVar += len(list_cVCF)

        OutFile.write('%s\t%s\n' % (sPatID, len(list_cVCF)))
        if sPatID not in dict_sVarDist:
            dict_sVarDist[sPatID] = {sLocale:0 for sLocale in list_sLocale}


        dict_sSNVType[sPatID] = [len(list_cVCF), len([cVCF for cVCF in list_cVCF if cVCF.sSNVType == 'nonsynonymous_SNV'])]

        dict_sLocaleKey = {'exonic':'CDS', '5UTR':'UTRs', '3UTR':'UTRs',
                           'intergenic':'intergenic', 'intronic':'intronic', 'other':'intergenic'}

        for cVCF in list_cVCF:
            #Split the Info column by ';' then key:value by '='
            dict_sVCFInfo = dict(sInfo.split('=') for sInfo in cVCF.sInfo.split(';')[:-1] if len(sInfo.split('=')) == 2)

            sKey          = dict_sGENIC_alt[cVCF.sLocale]
            sKey2         = cVCF.sGeneSym.split(',')[0] if ',' in cVCF.sGeneSym else cVCF.sGeneSym
            if sKey2 not in dict_sGeneSym:
                dict_sGeneSym[sKey2] = []
            dict_sGeneSym[sKey2].append(dict_sLocaleKey[sKey])
            dict_sGeneSym[sKey2] = sorted(list(set(dict_sGeneSym[sKey2])))

            if sKey not in dict_sVarDist[sPatID]:
                dict_sVarDist[sPatID][sKey] = 0
            dict_sVarDist[sPatID][sKey] += 1

            if sKey not in dict_sLocale:
                dict_sLocale[sKey] = 0
            dict_sLocale[sKey] += 1

            #list_tLOD.append(float(dict_sVCFInfo['TLOD']))
        #loop END: cVCF
    #loop END: sPatID
    OutFile.close()
    print(len(dict_sGeneSym))
    nCnt          = 0
    dict_sVarDist = {}
    for sGeneSym in dict_sRefSeq:

        try: dict_sGeneSym[sGeneSym]
        except KeyError:
            nCnt += 1
            continue

        sKey = '-'.join(dict_sGeneSym[sGeneSym])

        if sKey not in dict_sVarDist:
            dict_sVarDist[sKey] = []
        dict_sVarDist[sKey].append(sGeneSym)
    #loop END: sGeneSym
    print('No RefSeq', nCnt)
    for sKey in dict_sVarDist:
        print(sKey, len(dict_sVarDist[sKey]))

    sys.exit()



    OutFile = open('%s/NonsynVarCnt%s.txt' % (sWorkDir, '_cap' if bFiltered else ''), 'w')
    print(len(dict_sSNVType))
    for sPatID in dict_sSNVType:
        #print(sPatID, dict_sSNVType[sPatID][0], dict_sSNVType[sPatID][1])
        OutFile.write('%s\t%s\t%s\n' % (sPatID, dict_sSNVType[sPatID][0], dict_sSNVType[sPatID][1]))
    OutFile.close()


    nTotalVar2 = 0
    for sLocale in dict_sLocale:
        nTotalVar2 += dict_sLocale[sLocale]
        print(sLocale, (dict_sLocale[sLocale]))
        #print(sLocale, (dict_sLocale[sLocale] / dict_sGenicSize[sLocale]) * 1000000)

    print('Total', len(list_tLOD))
    print('MinLODt', min(list_tLOD))
    print('MAXLODt', max(list_tLOD))

    # VarDist by Patient
    list_sLocale = ['exonic', '5UTR', '3UTR', 'intergenic', 'intronic']
    sOutFile     = '%s/VariantsDistPerPatients%s.txt' % (sWorkDir, '_cap' if bFiltered else '')
    OutFile2     = open(sOutFile, 'w')
    for sPatID in dict_sVarDist:
        try: list_nVarCnts = [dict_sVarDist[sPatID][sLocale] for sLocale in list_sLocale]
        except KeyError:
            print(sPatID, dict_sVarDist[sPatID])
            sys.exit()


        sOut          = '%s\t%s\n' % (sPatID, '\t'.join([str(nCnt) for nCnt in list_nVarCnts]))
        OutFile2.write(sOut)
    #loop END: sPatID
    OutFile2.close()
#def END: basic_stats_vcf


def output_vcfdata_for_clustering (sWorkDir, bFiltered=0):
    print(green('Prepare VCFs for clustering', 1))

    sVCFDir         = '%s/VCFs'   % sWorkDir

    list_sPatIDs    = get_patientIDs(sSTUDY, sSEQDATA)
    list_cVCF_all   = []

    for sPatID in list_sPatIDs:
        print('Processing', sPatID)
        if bFiltered:
            sVCFFile = '%s/%s.hg19_multianno.dat' % (sVCFDir, sPatID)
            InFile = open(sVCFFile, 'rb')
            list_cVCF = pickle.load(InFile)
            InFile.close()
        else:
            sVCFFile = '%s/%s.hg19_multianno.vcf' % (sVCFDir, sPatID)
            print(sVCFFile)
            list_cVCF = cVCF_parse_vcf_files(sVCFFile)
            # if END:

        for cVCF in list_cVCF:
            cVCF.nVarCnt = len(list_cVCF)
            cVCF.sPatID  = sPatID
            list_cVCF_all.append(cVCF)
        #loop END: cVCF
    #loop END: sPatID
    list_cVCF_all = sorted(list_cVCF_all, key=lambda c:(c.nChrID, c.nPos))

    sOutFile = '%s/FullVCFsForClustering%s.txt' % (sWorkDir, '_cap' if bFiltered else '')
    OutFile  = open(sOutFile, 'w')
    for cVCF in list_cVCF_all:
        sOut = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
               % (cVCF.sChrID, cVCF.nPos, '.', cVCF.sRefNuc, cVCF.sAltNuc, cVCF.nVarCnt, cVCF.sInfo, cVCF.sPatID)
        OutFile.write(sOut)
    #loop END: cVCF
    OutFile.close()
#def END: output_vcfdata_for_clustering


def qsub_bedtools_cluster (sWorkDir, nClustDist, sQueue, bTestRun, bFiltered=0):
    print(green('STEP 2: Cluster using BEDtools', 1))

    sJobName = 'Jinman.BedTools.Cluster.%sbps'          % nClustDist
    sLogDir  = '%s/log/%s/%s'                           % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sLogDir, exist_ok=True)

    sVCFFile = '%s/FullVCFsForClustering%s.txt'         % (sWorkDir, '_cap' if bFiltered else '')
    sOutFile = '%s/FullVCFs.clustered.%sbps%s.txt'      % (sWorkDir, nClustDist, '_cap' if bFiltered else '')

    sScript  = 'bedtools cluster -d %s -i %s > %s'      % (nClustDist, sVCFFile, sOutFile)

    if bTestRun: print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s'
                  % (sScript, sLogDir, sQueue, sJobName))
    #if END
#def END: qsub_bedtools_cluster

## endregion

## region miRNA targets

def process_mir_data (sWorkDir, sRefFile):
    print(green('Process miR data', 1))

    print('Opening Pickled', sRefFile)
    OutFile         = open(sRefFile, 'rb')
    dict_cRefSeq    = pickle.load(OutFile)
    OutFile.close()

    dict_c3UTRs     = c3UTRs_load_3UTRS(s3UTR_FILE[sGENOME])

    #Assign 3UTR genomic position info
    annotate_3UTRs (dict_c3UTRs, dict_cRefSeq)

    if fTOP_EXPR != 1:
        #Load TCGA FPKM File List
        sFPKMDir        = '%s/%s/FPKMfiles/UQ_norm_recode'    % (sBASE_DIR, sSTUDY)
        sFPKMFile       = '%s/%s/%s_rnaseq-UQ.txt'            % (sBASE_DIR, sSTUDY, sSTUDY)
        dict_sFPKMFiles = load_fpkm_list (sFPKMFile)
        dict_sExpr      = get_expr_data (sFPKMDir, dict_sFPKMFiles)
    else:
        dict_sExpr      = {}

    # Obtain top miRNAs target sites in target cells (HepG2 or K562)
    list_cMir       = load_top_mirnas (dict_sMIRFILE[sCELLS])
    dict_sTargets   = determine_site_type_cnt (dict_c3UTRs, list_cMir, dict_sExpr)
    dict_sTargetPos = get_targetpos_by_gene (dict_sTargets)

    sOutFile        = '%s/mir_targets_top%s%%.dat' % (sWorkDir, int(fTOP_EXPR*100))
    OutFile         = open(sOutFile, 'wb')
    pickle.dump(dict_sTargetPos, OutFile)
    OutFile.close()
#def END: process_mir_data


def annotate_3UTRs (dict_c3UTRs, dict_cRefSeq):

    dict_cRef = {dict_cRefSeq[sGeneSym].sNMID:dict_cRefSeq[sGeneSym] for sGeneSym in dict_cRefSeq }

    for sNMID in dict_c3UTRs:

        c3UTR                = dict_c3UTRs[sNMID]
        cRef                 = dict_cRef[sNMID]

        c3UTR.sChrID         = cRef.sChrID
        c3UTR.sStrand        = cRef.sStrand
        c3UTR.list_n5UTRS    = cRef.list_n5UTRS
        c3UTR.list_n5UTRE    = cRef.list_n5UTRE
        c3UTR.list_nORFS     = cRef.list_nORFS
        c3UTR.list_nORFE     = cRef.list_nORFE
        c3UTR.list_n3UTRS    = cRef.list_n3UTRS
        c3UTR.list_n3UTRE    = cRef.list_n3UTRE
    #loop END: sNMID
#def END: annotate_3UTRs


def load_fpkm_list (sInFile):

    InFile       = open(sInFile, 'r')
    dict_sOutput = {}
    for sReadLine in InFile:
        ''' File Format
        0 sBarcode        A10Q
        1 sFPKMFile_T     e621bfc8-6557-450e-8e4a-3f754294cb20.FPKM-UQ.txt.gz
        2 sFPKMFile_N     576fda32-4e0c-414d-b54c-259f2503fc11.FPKM-UQ.txt.gz
        3 sPatID          TCGA_LIHC_Mutect2_349
        '''

        list_sColumn = sReadLine.strip('\n').split('\t')
        sBarCode     = list_sColumn[0]
        sFPKMFile_T  = list_sColumn[1]
        sFPKMFile_N  = list_sColumn[2]
        sPatID       = list_sColumn[3]

        if sPatID not in dict_sOutput:
            dict_sOutput[sPatID] = ''
        dict_sOutput[sPatID] = {'T':sFPKMFile_T, 'N':sFPKMFile_N, 'Barcode': sBarCode}
    #loop END:

    #VS-Check
    if not dict_sOutput:
        sys.exit('Invalid Dict : load_fpkm_list : dict_sOutput size= %d' % len(dict_sOutput))

    return dict_sOutput
#def END: load_fpkm_list


def get_expr_data (sFPKMDir, dict_sFPKMFiles):
    dict_sGeneSym      = {}
    for sPatID in dict_sFPKMFiles:
        sFPKMFile   = '%s/%s.N.txt' % (sFPKMDir, sPatID)
        dict_fFPKM  = load_fpkm(sFPKMFile)

        for sGeneSym in dict_fFPKM:
            if sGeneSym not in dict_sGeneSym:
                dict_sGeneSym[sGeneSym] = []
            dict_sGeneSym[sGeneSym].append(dict_fFPKM[sGeneSym])
        #loop END: sGeneSym
    #loop END: sPatID

    dict_sOutput = {}
    for sGeneSym in dict_sGeneSym:
        if sGeneSym not in dict_sOutput:
            dict_sOutput[sGeneSym] = 0.0
        dict_sOutput[sGeneSym] = np.mean(dict_sGeneSym[sGeneSym]) if sum(dict_sGeneSym[sGeneSym]) != 0 else 0.0
    #loop END: sGeneSym

    #VS-Check
    if not dict_sOutput:
        sys.exit('Invalid Dict : assign_fpkm_values : dict_sOutput size= %d' % len(dict_sOutput))

    return dict_sOutput
#def END: assign_fpkm_values


def load_top_mirnas (sInFile):
    list_cOutput = []
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:
        list_sColumn    = sReadLine.strip('\n').split('\t')

        cMir            = cMiRNA()
        cMir.sMirName   = list_sColumn[0]
        cMir.sMirSeq    = list_sColumn[1]
        cMir.dict_sST   = determine_site_types(cMir.sMirSeq)

        list_cOutput.append(cMir)
    #loop END: sReadLine
    return list_cOutput
#def END: load_top_mirnas


def determine_site_type_cnt (dict_c3UTRs, list_cMir, dict_sExpr):
    bSingleSite    = False
    dict_c3UTRs    = {dict_c3UTRs[sNMID].sGeneSym : dict_c3UTRs[sNMID] for sNMID in dict_c3UTRs}

    #Get Top Expressed
    #list_sTopGenes = [[sGeneSym, dict_sExpr[sGeneSym]] for sGeneSym in dict_sExpr if dict_sExpr[sGeneSym] > 0]
    #list_sTopGenes = sorted(list_sTopGenes, key=lambda e:e[1], reverse=True)[:int(len(list_sTopGenes)*fTOP_EXPR)]

    dict_sGeneSym  = {}

    for sGeneSym in dict_c3UTRs:

        try: c3UTR = dict_c3UTRs[sGeneSym]
        except KeyError: continue

        #Filter 3UTR Length
        if c3UTR.n3UTRLen > nMAX_3UTR_SIZE: continue

        #c3UTR.fExpr = fFPKM

        if sGeneSym not in dict_sGeneSym:
            dict_sGeneSym[sGeneSym] = ''
        dict_sGeneSym[sGeneSym] = c3UTR
    #loop END: sNMID

    list_c3UTRs   = sorted(list(dict_sGeneSym.values()), key=lambda c:c.n3UTRLen, reverse=True)

    dict_sTargets = {}
    for cMir in list_cMir:

        for c3UTR in list_c3UTRs:

            c3UTR             = copy.deepcopy(c3UTR)
            c3UTR.dict_nSTCnt = {sST : 0  for sST in list_sAllST}
            c3UTR.dict_nSTPos = {sST : [] for sST in list_sAllST}
            c3UTR.sMirName    = cMir.sMirName
            c3UTR.sMirSeq     = cMir.sMirSeq
            sTempSeq          = c3UTR.s3UTRSeq

            for sST in list_sCST:
                sSiteSeq  = cMir.dict_sST[sST]
                sTempSeq, c3UTR.dict_nSTCnt[sST], list_nPos = count_site_types(sTempSeq, sSiteSeq)

                for nIndexS, nIndexE in list_nPos:
                    if c3UTR.sStrand == '+':
                        nStart = c3UTR.list_n3UTRS[0] + nIndexS + 1
                        nEnd   = c3UTR.list_n3UTRS[0] + nIndexE + 1
                        c3UTR.dict_nSTPos[sST].append([nStart, nEnd])

                    else:
                        nStart = c3UTR.list_n3UTRE[-1] - nIndexE
                        nEnd   = c3UTR.list_n3UTRE[-1] - nIndexS
                        c3UTR.dict_nSTPos[sST].append([nStart, nEnd])
                #loop END: nIndexS, nIndexE

            #loop END: sST

            nAllST    = sum([c3UTR.dict_nSTCnt[sST] for sST in list_sAllST if sST != 'NC'])
            nCSTCnt   = sum([c3UTR.dict_nSTCnt[sST] for sST in list_sCST])
            nNSTCnt   = sum([c3UTR.dict_nSTCnt[sST] for sST in list_sNST])
            nCDNSTCnt = sum([c3UTR.dict_nSTCnt[sST] for sST in list_sCDNST if sST != 'NC'])

            # Determine UTRS with CSTs with Single or Multi Condition
            for sST in list_sAllST:

                if bSingleSite == True:

                    if c3UTR.dict_nSTCnt[sST] == 1 and sum([c3UTR.dict_nSTCnt[sWO_ST] for sWO_ST in dict_WO_ST[sST]]) == 0:
                        c3UTR.bSTFlag = sST

                    elif c3UTR.dict_nSTCnt[sST] >= 1 and sum([c3UTR.dict_nSTCnt[sWO_ST] for sWO_ST in dict_WO_ST[sST]]) == 0:
                        c3UTR.bSTFlag = 'Multi'

                else:

                    if c3UTR.dict_nSTCnt[sST] >= 1 and sum([c3UTR.dict_nSTCnt[sWO_ST] for sWO_ST in dict_WO_ST[sST]]) == 0:
                        c3UTR.bSTFlag = sST
                #if END: bSingleSite
            #loop END: sST

            if nAllST == 0:                                      # No Site
                c3UTR.bSTFlag = 'NoSite'
                cMir.list_c3UTRs_NoSite.append(c3UTR)

            elif nCSTCnt > 0 and (nNSTCnt + nCDNSTCnt) == 0:     # CST
                cMir.list_c3UTRs_CST.append(c3UTR)

            elif nCSTCnt > 0 and (nNSTCnt + nCDNSTCnt) > 0:      # CST
                cMir.list_c3UTRs_CST.append(c3UTR)

            elif nNSTCnt > 0 and (nCSTCnt + nCDNSTCnt) == 0:     # NST
                cMir.list_c3UTRs_NST.append(c3UTR)

            elif nNSTCnt > 0 and nCSTCnt == 0 and nCDNSTCnt > 0: # NST
                cMir.list_c3UTRs_NST.append(c3UTR)

            elif nCDNSTCnt > 0 and (nCSTCnt + nNSTCnt) == 0:     # CDNST
                cMir.list_c3UTRs_CDNST.append(c3UTR)
            else:
                print('What is This??', c3UTR.sNMID, nAllST, nCSTCnt, nNSTCnt, nCDNSTCnt)
            #if END: nCSTCnt vs nNSTCnt
        #loop END: c3UTR

        cMir.list_c3UTRs_NoSite = [c3UTR for c3UTR in cMir.list_c3UTRs_NoSite if c3UTR.n3UTRLen >= nMIN_NO_SITE_SIZE]

        if cMir.sMirName not in dict_sTargets:
            dict_sTargets[cMir.sMirName] = ''
        dict_sTargets[cMir.sMirName] = cMir
        #loop END: cMir

    return dict_sTargets
#def END: determine_site_type_cnt


def get_targetpos_by_gene (dict_sTargets):

    dict_sTargetPos = {}
    for sMirName in dict_sTargets:
        cMir = dict_sTargets[sMirName]
        for c3UTR in cMir.list_c3UTRs_CST:

            if c3UTR.sGeneSym not in dict_sTargetPos:
                dict_sTargetPos[c3UTR.sGeneSym] = {}

            for sST in list_sCST:
                if not c3UTR.dict_nSTPos[sST]: continue

                for nStart, nEnd in c3UTR.dict_nSTPos[sST]:

                    sKey = '%s:%s-%s' % (c3UTR.sChrID, nStart, nEnd)

                    if sKey not in dict_sTargetPos[c3UTR.sGeneSym]:
                        dict_sTargetPos[c3UTR.sGeneSym][sKey] = []
                    dict_sTargetPos[c3UTR.sGeneSym][sKey].append([sMirName, sST])
                #loop END: nStart, nEnd
            #loop END: sST
        #loop END: c3UTR
    #loop END: sMirName

    return dict_sTargetPos
#def END: get_targetpos_by_gene

## endregion

## endregion

## region Capture filter

def capture_filter (nClustDist, nBufferSize, sTempScript, sQueue, bTestRun):
    print(blue('Capture filter', 1))
    sPrevDir       = '%s/01_preprocess_data'            % sBASE_DIR
    sWorkDir       = '%s/02_capture_filter/%s'          % (sBASE_DIR, sSTUDY)
    sRefDir        = '%s/02_capture_filter/ref_capture' % sBASE_DIR
    sClipDir       = '%s/clipdata_forTL'                % sWorkDir
    sVCFDir        = '%s/VCFs-%s'                       % (sWorkDir, sSEQDATA)
    sMirDir        = '%s/mir'                           % sWorkDir
    os.makedirs(sWorkDir, exist_ok=True)
    os.makedirs(sRefDir, exist_ok=True)
    os.makedirs(sClipDir, exist_ok=True)
    os.makedirs(sVCFDir, exist_ok=True)
    os.makedirs(sMirDir, exist_ok=True)

    #Filter RefSeq
    filter_refseq (sRefDir, sPrevDir)

    #Filter eCLIP data
    #filter_eclip (sClipDir, sPrevDir, sTempScript, sQueue, bTestRun)

    #Filter VCF
    filter_vcf (sVCFDir, sPrevDir, nClustDist, sTempScript, sQueue, bTestRun)

    #Filter miRNAs
    #filter_miRNA_targets (sMirDir, sPrevDir)
#def END: capture_filter

## region Filter RefSeq

def filter_refseq (sWorkDir, sPrevDir):
    print(green('Filter RefSeq data', 1))
    #Load Capture bedfile
    list_sCapture = parse_capture_bedfile(sCAPTUREFILE)

    sCaptureName  = sCAPTUREFILE.split('/')[-1].split('.')[0]
    sOutFile      = '%s/refseq_%s.dat' % (sWorkDir, sCaptureName)

    if not os.path.isfile(sOutFile):
        #Load RefSeq
        sRefFile       = '%s/refseq.dat' % sPrevDir
        print('Opening Pickled', sRefFile)
        OutFile        = open(sRefFile, 'rb')
        dict_cRefSeq   = pickle.load(OutFile)
        OutFile.close()

        cGenome        = cFasta(dict_sGENOME[sGENOME])

        dict_sCapGene  = list_sCapture[0]
        dict_sCapID    = list_sCapture[1]
        print(len(dict_sCapGene), len(dict_sCapID))

        nCnt           = 0
        nCnt1          = 0

        list_sFiltered = []
        for sGeneSym in dict_cRefSeq:
            cRef      = dict_cRefSeq[sGeneSym]
            dict_sPos = {}

            # GeneSym and NMID do not match on same chrID
            try: dict_sPos.update(dict_sCapGene[cRef.sGeneSym])
            except KeyError: dict_sPos.update({})

            try: dict_sPos.update(dict_sCapID[cRef.sNMID])
            except KeyError: dict_sPos.update({})

            if not dict_sPos:
                nCnt += 1
                continue # Genes not found in capture file

            bCovered = cRef_determine_capture_coverage (cRef, dict_sPos)
            if not bCovered:
                nCnt1 += 1
                continue # Genes found but not covered in capture file

            cRef_update_seqinfo (cRef, cGenome)
            list_sFiltered.append(cRef)
        #if END: bCapture
        print('Not Found %s\tNot Covered %s\t Captured Genes %s'  % (nCnt, nCnt1, len(list_sFiltered)))

        OutFile        = open(sOutFile, 'wb')
        pickle.dump(list_sFiltered, OutFile)
        OutFile.close()
    else:
        print('%s FOUND' % sOutFile)
#def END: filter_refseq

## endregion

## region Filter eCLIP

def filter_eclip (sWorkDir, sPrevDir, sTempScript, sQueue, bTestRun):
    print(green('Filter eCLIP data', 1))

    sRefFile       = '%s/refseq.dat'    % sPrevDir
    sPrevDir       = '%s/%s'            % (sPrevDir, sSTUDY)
    sPrevDir       = '%s/clipdir_forTL' % sPrevDir
    sClipBedDir    = '%s/sorted_bed'    % sECLIP_DIR
    list_cClipData   = get_clipdata_list(1)

    list_sTargetRBPs = pathway_survey (list_cClipData)
    list_cClipData   = [cClipData for cClipData in list_cClipData if cClipData.sAccID in list_sTargetRBPs]
    print(len(list_sTargetRBPs))
    print(len(list_cClipData))

    '''Filter eCLIP with Capture'''
    #qsub_filter_eclipdata (sRefFile, sPrevDir, list_cClipData, sClipBedDir, sWorkDir, sTempScript, sQueue, bTestRun)

    '''Basic Stats of Filtered eCLIP data'''
    get_peak_distribution (sWorkDir, list_cClipData, 1)
#def END: filter_eclip


def qsub_filter_eclipdata (sRefFile, sPrevDir, list_cClipData, sClipBedDir, sWorkDir, sTempScript, sQueue, bTestRun):
    sJobName        = 'Jinman.CaptureFilter.eCLIP.%s' % sCELLS
    sLogDir         = '%s/log/%s/%s'                  % (sBASE_DIR, sJobName, sTIME_STAMP)
    sOutDir         = '%s/peaks'                      % sWorkDir
    sOutbyChr       = '%s/peaks_bychr'                % sWorkDir
    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sOutbyChr, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    for cClipData in list_cClipData:
        sBedFile       = '%s/%s.bed'                           % (sClipBedDir, cClipData.sAccID)
        sOutFile       = '%s/%s.peakdata'                      % (sOutDir, cClipData.sAccID)
        sScript        = '%s parse_eclip_stdout %s %s %s %s %s '  \
                         % (sTempScript, cClipData.sAccID, sRefFile, sBedFile, sOutFile, sOutbyChr)
        if bTestRun: print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s'
                      % (sScript, sLogDir, sQueue, sJobName, cClipData.sAccID))
        #if END
    #loop END: cClipData
    sBedFile       = '%s/clipdata/%s_combined.sorted.merged.bed' % (sPrevDir, sCELLS)
    sOutFile       = '%s/%s_combined.sorted.merged.filter.dat'   % (sWorkDir, sCELLS)
    sScript        = '%s parse_eclip_stdout_alt %s %s %s ' % (sTempScript, sRefFile, sBedFile, sOutFile)
    if bTestRun: print(sScript)
    else:
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.Merged'
                  % (sScript, sLogDir, sQueue, sJobName))
    #if END:
#def END: qsub_filter_eclipdata


def parse_eclip_stdout (sAccID, sRefFile, sInFile, sOutFile, sOutbyChr):
    print('Opening Pickled', sRefFile)
    OutFile        = open(sRefFile, 'rb')
    dict_cRefSeq   = pickle.load(OutFile)
    OutFile.close()

    dict_cRef      = {}
    for sNMID in dict_cRefSeq:
        sKey = '%s,%s' % (dict_cRefSeq[sNMID].sChrID, dict_cRefSeq[sNMID].sStrand)
        if sKey not in dict_cRef:
            dict_cRef[sKey] = []
        dict_cRef[sKey].append(dict_cRefSeq[sNMID])
    #loop END: sNMID

    dict_sChrID    = {sChrID:[] for sChrID in list_sCHRIDs}
    sScript        = 'bedtools intersect -a %s -b %s -wo'  % (sInFile, sCAPTUREFILE)
    sStdOut        = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout
    list_sOutput   = []
    for sReadLine in sStdOut:
        #Standard out format
        '''
        0   chr1                        sChrID                  Peak Info
        1   880079                      nPeakStart
        2   880152                      nPeakEnd
        3   GRSF1_HepG2_rep02           sSample
        4   200                         nScore
        5   -                           sStrand
        6   2.74834767396973            fSignal
        7   2.63395079019976            fPvalue
        8   -1                          NOT USED
        9   -1                          NOT USED
        10  chr1                        sChrID_cap              Capture Info
        11  880057                      nCapStart
        12  880208                      nCapEnd
        13  refGene=....                sGeneInfo
        14  73                          nOverlap
        '''
        sReadLine            = str(sReadLine, 'UTF-8')
        list_sColumn         = sReadLine.strip('\n').split('\t')

        cPeak                = cPeakData()
        cPeak.sChrID         = list_sColumn[0]
        cPeak.nStartPos      = max(int(list_sColumn[1]), int(list_sColumn[11]))
        cPeak.nEndPos        = min(int(list_sColumn[2]), int(list_sColumn[12]))

        assert cPeak.nStartPos < cPeak.nEndPos

        cPeak.nPeakSize      = cPeak.nEndPos - cPeak.nStartPos
        cPeak.nCenterPos     = cPeak.nEndPos - int(cPeak.nPeakSize / 2)

        cPeak.sSample        = list_sColumn[3]
        cPeak.nScore         = int(list_sColumn[4])
        cPeak.sStrand        = list_sColumn[5]
        cPeak.fSignal        = list_sColumn[6]
        cPeak.fPvalue        = list_sColumn[7]
        cPeak.list_sGeneInfo = [sInfo.split('=') for sInfo in list_sColumn[13].split(';') if sInfo] if list_sColumn[13] != '.' else []
        cPeak.nCapOverlap    = int(list_sColumn[14])

        determine_peak_locale(cPeak, dict_cRef)

        dict_sChrID[cPeak.sChrID].append(cPeak)

        list_sOutput.append(cPeak)

    #loop END: sReadLine

    #VS-Check
    if not list_sOutput:
        sys.exit('Invalid List : parse_eclip_stdout : list_sOutput size= %d' % len(list_sOutput))

    print(sOutFile, len(list_sOutput))

    OutFile         = open(sOutFile, 'wb')
    pickle.dump(list_sOutput, OutFile)
    OutFile.close()

    for sChrID in list_sCHRIDs:

        print(sChrID, len(dict_sChrID[sChrID]))

        sOutChr = '%s/%s' % (sOutbyChr, sChrID)
        os.makedirs(sOutChr, exist_ok=True)

        OutFile         = open('%s/%s.%s.peakdata' % (sOutChr, sChrID, sAccID), 'wb')
        pickle.dump(dict_sChrID[sChrID], OutFile)
        OutFile.close()
    #loop END: sChrID
#def END :parse_eclip_stdout


def parse_eclip_stdout_alt (sRefFile, sInFile, sOutFile):
    print('Opening Pickled', sRefFile)
    OutFile        = open(sRefFile, 'rb')
    dict_cRefSeq   = pickle.load(OutFile)
    OutFile.close()

    dict_cRef      = {}
    for sNMID in dict_cRefSeq:
        sKey = '%s,%s' % (dict_cRefSeq[sNMID].sChrID, dict_cRefSeq[sNMID].sStrand)
        if sKey not in dict_cRef:
            dict_cRef[sKey] = []
        dict_cRef[sKey].append(dict_cRefSeq[sNMID])
    #loop END: sNMID

    sScript        = 'bedtools intersect -a %s -b %s -wo'  % (sInFile, sCAPTUREFILE)
    sStdOut        = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout
    list_sOutput   = []
    for sReadLine in sStdOut:
        #Standard out format
        '''
        0   chr1                        sChrID                  Peak Info
        1   880079                      nPeakStart
        2   880152                      nPeakEnd
        3   GRSF1_HepG2_rep02           sSample
        4   200                         nScore
        5   -                           sStrand
        6   2.74834767396973            fSignal
        7   2.63395079019976            fPvalue
        8   -1                          NOT USED
        9   -1                          NOT USED
        10  chr1                        sChrID_cap              Capture Info
        11  880057                      nCapStart
        12  880208                      nCapEnd
        13  refGene=....                sGeneInfo
        14  73                          nOverlap
        '''
        sReadLine            = str(sReadLine, 'UTF-8')
        list_sColumn         = sReadLine.strip('\n').split('\t')

        cPeak                = cPeakData()
        cPeak.sChrID         = list_sColumn[0]
        cPeak.nStartPos      = max(int(list_sColumn[1]), int(list_sColumn[11]))
        cPeak.nEndPos        = min(int(list_sColumn[2]), int(list_sColumn[12]))

        assert cPeak.nStartPos < cPeak.nEndPos

        cPeak.nPeakSize      = cPeak.nEndPos - cPeak.nStartPos
        cPeak.nCenterPos     = cPeak.nEndPos - int(cPeak.nPeakSize / 2)

        cPeak.list_sSample   = list_sColumn[3].split(',')
        cPeak.list_nScore    = list_sColumn[4].split(',')
        cPeak.sStrand        = list_sColumn[5].split(',')[0]
        cPeak.list_fSignal   = list_sColumn[6].split(',')
        cPeak.list_fPvalue   = list_sColumn[7].split(',')
        cPeak.list_sGeneInfo = [sInfo.split('=') for sInfo in list_sColumn[13].split(';') if sInfo] if list_sColumn[13] != '.' else []
        cPeak.nCapOverlap    = int(list_sColumn[14])

        determine_peak_locale(cPeak, dict_cRef)

        list_sOutput.append(cPeak)
    #loop END: sReadLine

    #VS-Check
    if not list_sOutput:
        sys.exit('Invalid List : parse_eclip_stdout : list_sOutput size= %d' % len(list_sOutput))

    print(sOutFile, len(list_sOutput))

    OutFile         = open(sOutFile, 'wb')
    pickle.dump(list_sOutput, OutFile)
    OutFile.close()
#def END :parse_eclip_stdout

## endregion

## region Filter VCF

def filter_vcf (sWorkDir, sPrevDir, nClustDist, sTempScript, sQueue, bTestRun):
    print(green('Filter VCF data', 1))
    sPrevDir = '%s/%s' % (sPrevDir, sSTUDY)

    '''Filter VCFs with bedtools'''
    qsub_filter_vcf (sWorkDir, sPrevDir, sTempScript, sQueue, bTestRun)

    '''Basic Stats of Filtered VCFs'''
    #qsub_basic_stats_vcf           (sWorkDir, sTempScript, sQueue, bTestRun)  # Annotate VCFs
    #basic_stats_vcf                (sWorkDir, 1)

    '''Prepare VCFs for clustering'''
    #output_vcfdata_for_clustering   (sWorkDir, 1)

    '''STEP 2: Cluster using BEDtools'''
    #qsub_bedtools_cluster           (sWorkDir, nClustDist, sQueue, bTestRun, 1)

#def END: filter_vcf


def qsub_filter_vcf (sWorkDir, sPrevDir, sTempScript, sQueue, bTestRun):

    sVCFDir     = '%s/VCFs-%s/VCFs'             % (sPrevDir, sSEQDATA)
    sOutDir     = '%s/VCFs'                     % sWorkDir
    sPrevJob    = 'Jinman.ANNOVAR.VCF.%s.%s'    % (sSEQDATA, sSTUDY)
    sJobName    = 'Jinman.CaptureFilter.VCF.%s' % sSTUDY
    sLogDir     = '%s/log/%s/%s'                % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sOutDir, exist_ok=True)
    #os.makedirs(sChrIDDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    # Load Patient IDs
    list_sPatIDs = get_patientIDs(sSTUDY, sSEQDATA)

    for sPatID in list_sPatIDs:
        sInFile  = '%s/%s.hg19_multianno.vcf'    % (sVCFDir, sPatID)
        sZipFile = '%s/%s.hg19_multianno.vcf.gz' % (sOutDir, sPatID)

        if os.path.isfile(sZipFile): continue
        sScript  = 'cat %s | egrep PASS | awk \'{ if (\\$2 != \\\"PASS\\\") print }\' |  bgzip -c > %s;' % (sInFile, sZipFile)
        sScript += 'tabix -p vcf %s;'  % sZipFile

        if bTestRun:
            print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s -hold_jid %s.%s'
                      % (sScript, sLogDir, sQueue, sJobName, sPatID, sPrevJob, sPatID))
    #loop END: sPatID
#def END: qsub_filter_vcf


def parse_vcf_stdout (sPatID, sInFile, sOutFile, sOutDir):

    if sSEQDATA == 'WXS':
        sScript        = 'bedtools intersect -a %s -b %s -wo'  % (sInFile, sCAPTUREFILE)
        sStdOut        = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout
    else:
        sStdOut        = open(sInFile, 'r')

    dict_sChrID    = {}
    list_sOutput   = []
    for sReadLine in sStdOut:
        # File Format
        # Column Number:     | 0       | 1        | 2          | 3       | 4
        # Column Description:| sChrID  | nPos     | sDBSNP_ID  | sRefNuc | sAltNuc
        # Column Example:    | chr13   | 32906558 | rs79483201 | T       | A
        # Column Number:     | 5       | 6        | 7          | 8              | 9./..
        # Column Description:| fQual   | sFilter  | sInfo      | sFormat        | sSampleIDs
        # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to sFormat
        ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph
        if sSEQDATA == 'WXS': sReadLine           = str(sReadLine, 'UTF-8')
        list_sColumn        = sReadLine.strip('\n').split('\t')

        cVCF                = cVCFData()

        if sSEQDATA == 'WXS':
            cVCF.sChrID  = list_sColumn[0]
            cVCF.nChrID  = process_chromo_ID(list_sColumn[0])
        else:
            if list_sColumn[0] == 'MT': continue
            if list_sColumn[0].startswith('<GL00'): continue

            dict_sChrKey = {'X':'23', 'Y':'24'}
            cVCF.sChrID  = 'chr%s' % list_sColumn[0]

            if list_sColumn[0] in ['X', 'Y']:
                cVCF.nChrID  = int(dict_sChrKey[list_sColumn[0]])
            else:
                cVCF.nChrID  = int(list_sColumn[0])
        #if END:

        cVCF.nPos           = int(list_sColumn[1])
        cVCF.sDBSNP_ID      = list_sColumn[2]
        cVCF.sRefNuc        = list_sColumn[3]
        cVCF.sAltNuc        = list_sColumn[4]
        cVCF.fQual          = float(list_sColumn[5]) if list_sColumn[5] != '.' else list_sColumn[5]
        cVCF.sFilter        = list_sColumn[6]
        cVCF.sInfo          = list_sColumn[7]
        cVCF.sFormat        = list_sColumn[8]
        cVCF.list_sSamples  = list_sColumn[9:]

        dict_sVCFInfo       = dict(sInfo.split('=') for sInfo in cVCF.sInfo.split(';')[:-1] if len(sInfo.split('=')) == 2)
        cVCF.sLocale        = dict_sGENIC_alt[dict_sVCFInfo['Func.refGene'].split('\\x3b')[0]]
        cVCF.sGeneSym       = dict_sVCFInfo['Gene.refGene'].split('\\x3b')[0].upper()
        cVCF.sSNVType       = dict_sVCFInfo['ExonicFunc.refGene'].split('\\x3b')[0]

        sKey                = cVCF.sChrID
        if sKey not in dict_sChrID:
            dict_sChrID[sKey] = []
        dict_sChrID[sKey].append(cVCF)

        if sSEQDATA == 'WXS':
            if cVCF.sFilter != 'PASS': continue
        else:
            if 'PASS' not in sReadLine:continue

        list_sOutput.append(cVCF)
    #loop END: sReadLine

    if not list_sOutput:
        sys.exit('Invalid List : parse_eclip_stdout : list_sOutput size= %d' % len(list_sOutput))

    print(sOutFile, len(list_sOutput))
    OutFile         = open(sOutFile, 'wb')
    pickle.dump(list_sOutput, OutFile)
    OutFile.close()

    for sChrID in dict_sChrID:

        sChrDir = '%s/%s' % (sOutDir, sChrID)
        os.makedirs(sChrDir, exist_ok=True)

        print(sChrID, len(dict_sChrID[sChrID]))

        sOutFile = '%s/%s.%s.hg19_multianno.dat' % (sChrDir, sPatID, sChrID)
        OutFile  = open(sOutFile, 'wb')
        pickle.dump(dict_sChrID[sChrID], OutFile)
        OutFile.close()
    #loop END: sChrID
#def END: parse_vcf_stdout

## endregion

## region Filter miRNAs target sites
def filter_miRNA_targets (sWorkDir, sPrevDir):

    print(green('Filter miR data', 1))
    sPrevDir        = '%s/%s' % (sPrevDir, sSTUDY)

    #Load miR target file
    sInFile         = '%s/mir/mir_targets_top%s%%.dat' % (sPrevDir, int(fTOP_EXPR*100))
    InFile          = open(sInFile, 'rb')
    dict_sTargetPos = pickle.load(InFile)
    print('dict_sTargetPos', len(dict_sTargetPos))
    InFile.close()

    #Load Capture bedfile
    list_sCapture  = parse_capture_bedfile(sCAPTUREFILE)
    dict_sCapGene  = list_sCapture[0]
    dict_sCapID    = list_sCapture[1]

    dict_c3UTRs    = c3UTRs_load_3UTRS(s3UTR_FILE[sGENOME])
    dict_c3UTRs    = {dict_c3UTRs[sNMID].sGeneSym : dict_c3UTRs[sNMID] for sNMID in dict_c3UTRs}

    dict_sFiltered = {}

    nCnt   = 0
    nCnt1  = 0
    for sGeneSym in dict_sTargetPos:

        c3UTR     = dict_c3UTRs[sGeneSym]

        dict_sPos = {}
        # GeneSym and NMID do not match on same chrID
        try: dict_sPos.update(dict_sCapGene[sGeneSym])
        except KeyError: dict_sPos.update({})

        try: dict_sPos.update(dict_sCapID[c3UTR.sNMID])
        except KeyError: dict_sPos.update({})

        if not dict_sPos:
            nCnt += 1
            continue # Genes not found in capture file

        dict_nNewPos = filter_capture_mir (dict_sTargetPos[sGeneSym], dict_sPos)

        if not dict_nNewPos:
            nCnt1 += 1
            continue # Genes not covered in capture file

        for sChrID in dict_nNewPos:
            if sChrID not in dict_sFiltered:
                dict_sFiltered[sChrID] = {}

            if sGeneSym not in dict_sFiltered[sChrID]:
                dict_sFiltered[sChrID][sGeneSym] = {}

            dict_sFiltered[sChrID][sGeneSym] = dict_nNewPos[sChrID]
        #loop END: sChrID
    #loop END: sGeneSym
    print('Not Found %s\tNot Covered %s\t Captured Genes %s' % (nCnt, nCnt1, len(dict_sFiltered)))

    sOutFile = '%s/mir_targets_top%s%%_cap.dat' % (sWorkDir, int(fTOP_EXPR*100))
    OutFile  = open(sOutFile, 'wb')
    pickle.dump(dict_sFiltered, OutFile)
    OutFile.close()
#def END: filter_miRNA_targets


def filter_capture_mir (dict_sTargetPos, dict_sPos):

    dict_sFiltered = {}
    for sPosKey in dict_sTargetPos:

        sChrID       = sPosKey.split(':')[0]
        nStart, nEnd = [int(sPos) for sPos in sPosKey.split(':')[1].split('-')]

        list_sCapPos = dict_sPos[sChrID] # VS-Check

        assert len(list(set([sST for sMirName, sST in dict_sTargetPos[sPosKey]]))) == 1
        sTargetST    = list(set([sST for sMirName, sST in dict_sTargetPos[sPosKey]]))[0]

        if sTargetST not in list_s78MER: continue

        for sCapKey in list_sCapPos:

            nCapStart, nCapEnd = [int(sPos) for sPos in sCapKey.split(',')]

            if nCapStart > nEnd: continue
            if nCapEnd < nStart: continue

            nNewStart = max(nStart, nCapStart)
            nNewEnd   = min(nEnd, nCapEnd)

            sKey = '%s-%s' % (nNewStart, nNewEnd)

            if sChrID not in dict_sFiltered:
                dict_sFiltered[sChrID] = {}

            if sKey not in dict_sFiltered[sChrID]:
                dict_sFiltered[sChrID][sKey] = []
            dict_sFiltered[sChrID][sKey] = dict_sTargetPos[sPosKey]
        #loop END: sCapKey
    #loop END: sPosKey

    return dict_sFiltered
#def END: filter_capture_mir
## endregion

## endregion

## region Determine Enrichment
def determine_enrich (nVarWindow, nClustDist, nBufferSize, fRecurrency, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun):
    print(blue('Determine enrichment', 1))

    nChrWin      = 1000000
    sPrepDir     = '%s/01_preprocess_data'     % sBASE_DIR
    sCapDir      = '%s/02_capture_filter'      % sBASE_DIR
    sWorkDir     = '%s/03_determine_enrich'    % sBASE_DIR

    sMirDir      = '%s/mir'                    % sPrepDir
    os.makedirs(sMirDir,   exist_ok=True)

    #Load RefSeq data
    list_cRefSeq = load_refseq_data (sCapDir if sSEQDATA == 'WXS' else sPrepDir)
    print('list_cRefSeq', len(list_cRefSeq))

    #Acquire top expressed genes
    #dict_sGeneSym = acquire_top_expressed_genes (list_cRefSeq)
    #print('dict_sGeneSym', len(dict_sGeneSym))

    if sSEQDATA == 'WGS':
        #for WGS, use multiplexed qsub functions

        ### Collapse eCLIP data ###
        sClipDir = '%s/clipdata' % sWorkDir
        os.makedirs(sClipDir, exist_ok=True)
        qsub_collapse_eclip_data (sClipDir, sPrepDir, bStrictPeaks, bAllRBPs, sTempScript, '24_730.q', bTestRun)
        ###########################

        ### Load Clustered VCF Data ###
        sVCFDir   = '%s/%s/VCFs-%s'   % (sWorkDir, sSTUDY, sSEQDATA)
        os.makedirs(sVCFDir, exist_ok=True)
        #qsub_collapse_vcf_by_chrwin (sVCFDir, sPrepDir, nChrWin, sTempScript, sQueue, bTestRun)
        ###############################

        ### Assign variant location ###
        sStudyDir = '%s/%s'           % (sWorkDir, sSTUDY)
        os.makedirs(sStudyDir, exist_ok=True)
        #qsub_assign_variants         (sStudyDir, sClipDir, sMirDir, nBufferSize, nChrWin, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun)
        #qsub_assign_variants_byRBP   (sStudyDir, sClipDir, sMirDir, nBufferSize, nChrWin, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun)
        #qsub_combine_variants_byRBP  (sStudyDir, nChrWin, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun)
        ###############################

        ### Determine Region by RBPs ###
        #combined_rbpsize            (sClipDir, bStrictPeaks, bAllRBPs)
        ################################

        sStudyDir = '%s/%s'           % (sWorkDir, sSTUDY)
        os.makedirs(sStudyDir, exist_ok=True)

        ### Determine Genic Region Sizes ###
        #calculate_regionsize (sStudyDir, sCapDir, sPrepDir, list_cRefSeq, bStrictPeaks, bAllRBPs)
        ####################################

        ### Output Normalized Enrichment by Chr or RBP by RBP ###
        bPatCnt = 1
        #normalize_enrich_bychr (sStudyDir, nChrWin, fRecurrency, bPatCnt, bStrictPeaks, bAllRBPs)
        normalize_enrich_byrbp_v1 (sStudyDir, sClipDir, fRecurrency, bPatCnt, bStrictPeaks, bAllRBPs)
        #qsub_normalize_enrich_byrbp  (sStudyDir, sClipDir, fRecurrency, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun)
        #########################################################
        pass
    elif sSEQDATA == 'WXS':
        # Load Clustered VCF Data
        # dict_cVCF    = collapse_vcf_by_pos (sVCFDir, sCapDir, nVarWindow, 1)
        # print('dict_cVCF', len(dict_cVCF))

        # dict_nSizeRef = calculate_regionsize (sStudyDir, sCapDir, sPrepDir, list_cRefSeq, bStrictPeaks)
        pass
    #if END:
#def END: determine_enrich
'''
def determine_enrich_old (nVarWindow, nClustDist, nBufferSize, fRecurrency, bStrictPeaks, sTempScript, sQueue, bTestRun):
    print(blue('Determine enrichment', 1))

    nChrWin      = 1000000
    sPrev2Dir    = '%s/01_preprocess_data'     % sBASE_DIR
    sPrevDir     = '%s/02_capture_filter'      % sBASE_DIR
    sWorkDir     = '%s/03_determine_enrich'    % sBASE_DIR
    sStudyDir    = '%s/%s'                     % (sWorkDir, sSTUDY)
    sClipDir     = '%s/clipdata'               % sWorkDir
    sVCFDir      = '%s/VCFs-%s'                % (sStudyDir, sSEQDATA)
    sMirDir      = '%s/mir'                    % sPrev2Dir
    os.makedirs(sStudyDir, exist_ok=True)
    os.makedirs(sClipDir,  exist_ok=True)
    os.makedirs(sVCFDir,   exist_ok=True)
    os.makedirs(sMirDir,   exist_ok=True)

    #Load RefSeq data
    list_cRefSeq = load_refseq_data (sPrevDir if sSEQDATA == 'WXS' else sPrev2Dir)
    print('list_cRefSeq', len(list_cRefSeq))

    #Acquire top expressed genes
    dict_sGeneSym = acquire_top_expressed_genes (list_cRefSeq)
    print('dict_sGeneSym', len(dict_sGeneSym))

    if sSEQDATA == 'WXS':

        #Load Clustered VCF Data
        #dict_cVCF    = collapse_vcf_by_pos (sVCFDir, sPrevDir, nVarWindow, 1)
        #print('dict_cVCF', len(dict_cVCF))

        #dict_nSizeRef = calculate_regionsize (sStudyDir, sPrevDir, sPrev2Dir, list_cRefSeq, bStrictPeaks)

        pass

    elif sSEQDATA == 'WGS':
        #for WGS, use multiplexed qsub functions

        #Collapse eCLIP data
        qsub_collapse_eclip_data     (sClipDir, sPrev2Dir, bStrictPeaks, sTempScript, '24_730.q', bTestRun)

        #Load Clustered VCF Data
        #qsub_collapse_vcf_by_chrwin  (sVCFDir, sPrevDir, nChrWin, sTempScript, sQueue, bTestRun)

        #Assign variant location
        #qsub_assign_variants         (sStudyDir, sClipDir, sMirDir, nBufferSize, nChrWin, bStrictPeaks, sTempScript, sQueue, bTestRun)
        #qsub_assign_variants_byRBP   (sStudyDir, sClipDir, sMirDir, nBufferSize, nChrWin, bStrictPeaks, sTempScript, sQueue, bTestRun)
        #qsub_combine_variants_byRBP  (sStudyDir, nChrWin, bStrictPeaks, sTempScript, sQueue, bTestRun)
        #combined_rbpsize             (sClipDir, bStrictPeaks)

        #dict_nSizeRef = calculate_regionsize (sStudyDir, sPrevDir, sPrev2Dir, list_cRefSeq, bStrictPeaks)

        #normalize_enrich_bychr     (sStudyDir, nChrWin, fRecurrency, bStrictPeaks)
        #normalize_enrich_byrbp      (sStudyDir, sClipDir, dict_nSizeRef, fRecurrency, bStrictPeaks)
        pass
    #if END:

    sys.exit()
    #Load miRNA data
    dict_sMirs   = load_mirna_data (sPrevDir if sSEQDATA == 'WXS' else sPrev2Dir)
    print('dict_sMirs', len(dict_sMirs))

    #Load eCLIP data
    dict_cPeak, dict_sRBP = load_eclip_data (sClipDir, sPrevDir if sSEQDATA == 'WXS' else sPrev2Dir)
    print('dict_cPeak', len(dict_cPeak))
    print('dict_sRBP', len(dict_sRBP))

    #Acquire top expressed genes
    dict_sGeneSym = acquire_top_expressed_genes (list_cRefSeq)
    print('dict_sGeneSym', len(dict_sGeneSym))

    #Assign variant location
    dict_nAssignedVars = assign_variants (sWorkDir, dict_cVCF, dict_sMirs, dict_cPeak, dict_sGeneSym, nBufferSize, nVarWindow)
    print('dict_nAssignedVars', len(dict_nAssignedVars))

    #Calculate RBP and non-RBP region sizes
    dict_nSizeRef = calculate_regionsize (sWorkDir, sPrevDir, sPrev2Dir, list_cRefSeq)
    #calculate_regionsize_by_rbp (sWorkDir, dict_sRBP,  dict_sGeneSym, nBufferSize, nVarWindow)

    #Determined enrichment
    normalize_enrich (dict_cVCF, dict_nAssignedVars, fRecurrency)
    ## Survey by Individual RBPs
    #assign_variants_by_rbp      (sWorkDir, dict_cVCF, dict_sMirs, dict_sRBP, dict_sGeneSym, nBufferSize, nVarWindow)
    #output_enrich_by_rbp        (sWorkDir, dict_cVCF, dict_sRBP, dict_nAssignedVars, dict_nSizeRef, fRecurrency, nBufferSize, nVarWindow)
#def END: determine_enrich
'''

def qsub_collapse_vcf_by_chrwin (sWorkDir, sPrevDir, nChrWin, sTempScript, sQueue, bTestRun):
    sVCFDir     = '%s/%s/VCFs-%s/VCFs'          % (sPrevDir, sSTUDY, sSEQDATA)
    sChrDir     = '%s/VCFs_bychrwin'            % sWorkDir
    sPrevJob    = 'Jinman.CaptureFilter.VCF.%s' % sSTUDY
    sJobName    = 'Jinman.CollapseByPos.VCF.%s' % sSTUDY
    sLogDir     = '%s/log/%s/%s'                % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sChrDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    for sChrID in list_sCHRIDs:

        list_nWindow    = set_chr_range(sChrID, nChrWin)

        for i, sViewRange in enumerate(list_nWindow):

            sInterval   = '%s.%s' % (i, sViewRange)
            sOutFile    = '%s/vcfs_%s_%s.%s.dat' % (sChrDir, sSEQDATA, sChrID, sViewRange)

            sScript     = '%s collapse_vcf_by_pos_tabix %s %s %s %s'        \
                           % (sTempScript, sChrID, sViewRange, sVCFDir, sOutFile)

            if bTestRun:
                print(sScript)
            else:
                os.makedirs(sLogDir, exist_ok=True)
                os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s.%s.%s -hold_jid %s.*'
                          % (sScript, sLogDir, sQueue, sJobName, sChrID, sSEQDATA, sInterval, sPrevJob))

        #loop END: i, sViewRange
    #loop END: sChrID
#def END: qsub_collapse_vcf_by_chrwin


def qsub_collapse_vcf_by_chr (sWorkDir, sPrevDir, sTempScript, sQueue, bTestRun):

    sVCFDir     = '%s/VCFs-%s/VCFs'       % (sPrevDir, sSEQDATA)
    sChrDir     = '%s/VCFs_bychr'         % sWorkDir
    sJobName    = 'Jinman.CollapseByPos.VCF'
    sLogDir     = '%s/log/%s/%s'          % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sChrDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    for sChrID in list_sCHRIDs:

        sInDir      = '%s/%s' % (sVCFDir, sChrID)
        sOutFile    = '%s/vcfs_%s_%s.dat' % (sChrDir, sSEQDATA, sChrID)

        if os.path.isfile(sOutFile): continue

        sScript  = '%s collapse_vcf_by_pos_bychr %s %s %s '        \
                    % (sTempScript, sChrID, sInDir, sOutFile)

        if bTestRun:
            print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s.%s'
                      % (sScript, sLogDir, sQueue, sJobName, sChrID, sSEQDATA))
    #loop END: sChrID
#def END: qsub_collapse_vcf_by_chr


def collapse_vcf_by_pos_tabix (sChrID, sViewRange, sInDir, sOutFile):
    print(green('Collapse VCFs by Position by Chr', 1))

    list_sPatIDs       = get_patientIDs(sSTUDY, sSEQDATA)
    nWinStart, nWinEnd = sViewRange.split('-')
    list_cVCF_all      = []

    for sPatID in list_sPatIDs:
        sInFile    = '%s/%s.hg19_multianno.vcf.gz' % (sInDir, sPatID)
        sScript    = 'tabix %s %s:%s-%s'           % (sInFile, sChrID[-1], nWinStart, nWinEnd)
        sStdOut    = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout
        list_cVCF  = parse_vcf_stdout2(sStdOut)

        for cVCF in list_cVCF:
            cVCF.nVarCnt = len(list_cVCF)
            cVCF.sPatID  = sPatID
            list_cVCF_all.append(cVCF)
        #loop END: cVCF
    #loop END: sPatID

    list_cVCF_all = sorted(list_cVCF_all, key=lambda c:(c.nChrID, c.nPos))

    print(sOutFile, len(list_cVCF_all))

    ## Group by cluster and patient IDs
    dict_nPos = {}
    for cVCF in list_cVCF_all:
        if cVCF.nPos not in dict_nPos:
            dict_nPos[cVCF.nPos] = []

        dict_nPos[cVCF.nPos].append(cVCF)
    #loop END: cVCF

    OutFile  = open(sOutFile, 'wb')
    pickle.dump(dict_nPos, OutFile)
    OutFile.close()
#def END: collapse_vcf_by_pos_tabix


def parse_vcf_stdout2 (sStdOut):

    list_sOutput   = []
    for sReadLine in sStdOut:
        # File Format
        # Column Number:     | 0       | 1        | 2          | 3       | 4
        # Column Description:| sChrID  | nPos     | sDBSNP_ID  | sRefNuc | sAltNuc
        # Column Example:    | chr13   | 32906558 | rs79483201 | T       | A
        # Column Number:     | 5       | 6        | 7          | 8              | 9./..
        # Column Description:| fQual   | sFilter  | sInfo      | sFormat        | sSampleIDs
        # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to sFormat
        ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph
        sReadLine           = str(sReadLine, 'UTF-8')
        list_sColumn        = sReadLine.strip('\n').split('\t')

        cVCF                = cVCFData()
        if list_sColumn[0] == 'MT': continue
        if list_sColumn[0].startswith('<GL00'): continue

        dict_sChrKey = {'X':'23', 'Y':'24'}
        cVCF.sChrID  = 'chr%s' % list_sColumn[0]

        if list_sColumn[0] in ['X', 'Y']:
            cVCF.nChrID  = int(dict_sChrKey[list_sColumn[0]])
        else:
            cVCF.nChrID  = int(list_sColumn[0])
        cVCF.nPos           = int(list_sColumn[1])
        cVCF.sDBSNP_ID      = list_sColumn[2]
        cVCF.sRefNuc        = list_sColumn[3]
        cVCF.sAltNuc        = list_sColumn[4]
        cVCF.fQual          = float(list_sColumn[5]) if list_sColumn[5] != '.' else list_sColumn[5]
        cVCF.sFilter        = list_sColumn[6]
        cVCF.sInfo          = list_sColumn[7]
        cVCF.sFormat        = list_sColumn[8]
        cVCF.list_sSamples  = list_sColumn[9:]

        dict_sVCFInfo       = dict(sInfo.split('=') for sInfo in cVCF.sInfo.split(';')[:-1] if len(sInfo.split('=')) == 2)
        cVCF.sLocale        = dict_sGENIC_alt[dict_sVCFInfo['Func.refGene'].split('\\x3b')[0]]
        cVCF.sGeneSym       = dict_sVCFInfo['Gene.refGene'].split('\\x3b')[0].upper()
        cVCF.sSNVType       = dict_sVCFInfo['ExonicFunc.refGene'].split('\\x3b')[0]

        list_sOutput.append(cVCF)
    #loop END: sReadLine

    return list_sOutput
#def END: parse_vcf_stdout


def collapse_vcf_by_pos_bychr (sChrID, sInDir, sOutFile):
    print(green('Collapse VCFs by Position by Chr', 1))

    list_sPatIDs    = get_patientIDs(sSTUDY, sSEQDATA)
    list_cVCF_all   = []
    for sPatID in list_sPatIDs:
        sInFile     = '%s/%s.%s.hg19_multianno.dat' % (sInDir, sPatID, sChrID)
        InFile      = open(sInFile, 'rb')
        list_cVCF   = pickle.load(InFile)
        InFile.close()

        for cVCF in list_cVCF:
            cVCF.nVarCnt = len(list_cVCF)
            cVCF.sPatID  = sPatID
            list_cVCF_all.append(cVCF)
        #loop END: cVCF
    #loop END: sPatID
    list_cVCF_all = sorted(list_cVCF_all, key=lambda c:(c.nChrID, c.nPos))

    print(sChrID, len(list_cVCF_all))

    ## Group by cluster and patient IDs
    dict_nPos = {}
    for cVCF in list_cVCF_all:
        if cVCF.nPos not in dict_nPos:
            dict_nPos[cVCF.nPos] = []

        dict_nPos[cVCF.nPos].append(cVCF)
    #loop END: cVCF
    OutFile  = open(sOutFile, 'wb')
    pickle.dump(dict_nPos, OutFile)
    OutFile.close()
#def END: collapse_vcf_by_pos_bychr


def collapse_vcf_by_pos (sWorkDir, sPrevDir, nVarWindow, bFiltered=0):
    print(green('Collapse VCFs by Position', 1))

    sVCFFile     = '%s/vcfs_%s.dat'       % (sWorkDir, sSEQDATA)
    sVCFFile_win = '%s/vcfs_%s_win%s.dat' % (sWorkDir, sSEQDATA, nVarWindow)

    if not os.path.isfile(sVCFFile):

        print('Making Pickles', sVCFFile)

        sInVCFDir       = '%s/VCFs-%s/VCFs_bychr'   % (sPrevDir, sSEQDATA)
        list_sPatIDs    = get_patientIDs(sSTUDY, sSEQDATA)
        list_cVCF_all   = []

        for sPatID in list_sPatIDs:
            if bFiltered:
                sInFile = '%s/%s.hg19_multianno.dat' % (sInVCFDir, sPatID)
                InFile  = open(sInFile, 'rb')
                list_cVCF = pickle.load(InFile)
                InFile.close()
            else:
                sInFile = '%s/%s.hg19_multianno.vcf' % (sInVCFDir, sPatID)
                list_cVCF = cVCF_parse_vcf_files(sInFile)
            # if END:

            for cVCF in list_cVCF:
                cVCF.nVarCnt = len(list_cVCF)
                cVCF.sPatID  = sPatID

                if cVCF.sFilter != 'PASS': continue

                list_cVCF_all.append(cVCF)
            #loop END: cVCF
        #loop END: sPatID
        list_cVCF_all = sorted(list_cVCF_all, key=lambda c:(c.nChrID, c.nPos))

        ## Group by cluster and patient IDs
        dict_nPos = {}
        for cVCF in list_cVCF_all:

            if cVCF.sChrID not in dict_nPos:
                dict_nPos[cVCF.sChrID] = {}

            if cVCF.nPos not in dict_nPos[cVCF.sChrID]:
                dict_nPos[cVCF.sChrID][cVCF.nPos] = []

            dict_nPos[cVCF.sChrID][cVCF.nPos].append(cVCF)
        #loop END: cVCF

        OutFile = open(sVCFFile, 'wb')
        pickle.dump(dict_nPos, OutFile)
        OutFile.close()

        return dict_nPos

    else:
        InFile      = open(sVCFFile, 'rb')
        dict_nPos   = pickle.load(InFile)
        InFile.close()

        if nVarWindow:
            dict_sOutput = {}
            for sChrID in dict_nPos:

                if sChrID not in dict_sOutput:
                    dict_sOutput[sChrID] = {}

                #Group by window
                for nPosID, lCnt in itertools.groupby(dict_nPos[sChrID], key=lambda n: n // nVarWindow):

                    list_nPos = list(lCnt)
                    nKey      = int(nPosID * nVarWindow + int(nVarWindow / 2))
                    if nKey not in dict_sOutput[sChrID]:
                        dict_sOutput[sChrID][nKey] = []

                    for nVarPos in list_nPos:
                        dict_sOutput[sChrID][nKey] += dict_nPos[sChrID][nVarPos]
                    #loop END: nVarPos
                # loop END: nGroup, lCnt
            # loop END: sChrID

            OutFile = open(sVCFFile_win, 'wb')
            pickle.dump(dict_sOutput, OutFile)
            OutFile.close()

            dict_sPositions = dict_sOutput
            list_nPosCnt = []
            dict_sHisto  = {}
            for sChrID in dict_sPositions:
                for nPos in dict_sPositions[sChrID]:

                    list_nPosCnt.append(nPos)

                    sKey = len(dict_sPositions[sChrID][nPos])

                    if sKey not in dict_sHisto:
                        dict_sHisto[sKey] = 0
                    dict_sHisto[sKey] += 1
                #loop END: nPos
            #loop END: sChrID

            list_sKey = sorted(list(dict_sHisto.keys()), reverse=True)
            print(len(list_nPosCnt))
            print(max(list_sKey))
            print(min(list_sKey))
            for sKey in list_sKey:
                fRecurrency = '%0.2f%%' % ((sKey / nSAMPLE_CNT) * 100)
                fPercentage = '%0.2f%%' % (dict_sHisto[sKey]/len(list_nPosCnt)*100)
                print(sKey, fRecurrency, dict_sHisto[sKey], fPercentage)
            #loop END: sKey
            return dict_sOutput
        else:
            return dict_nPos
#def END: collapse_vcf_by_pos


def qsub_collapse_eclip_data (sClipDir, sPrepDir, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun):
    # Proofread 6/22/17
    sRBPGroup   = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall   = 'strict' if bStrictPeaks else 'allpeaks'
    sPeakDir    = '%s/clipdata/%s/%s-peaks_bychr'          % (sPrepDir, sPeakCall, sRBPGroup)
    sChrDir     = '%s/%s/%s-RBPs_bychr'                    % (sClipDir, sPeakCall, sRBPGroup)
    sSizeDir    = '%s/%s/%s-RBPsizes'                      % (sClipDir, sPeakCall, sRBPGroup)
    sJobName    = 'Jinman.CollapseByPos.eCLIP.%s.%s.%s.%s' % (sSTUDY, sSEQDATA, sPeakCall, sRBPGroup)
    sLogDir     = '%s/log/%s/%s'                           % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sChrDir, exist_ok=True)
    os.makedirs(sSizeDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    for sChrID in list_sCHRIDs:
        sInDir   = '%s/%s'              % (sPeakDir, sChrID)
        sOutFile = '%s/peaks_%s_%s.dat' % (sChrDir, sSEQDATA, sChrID)

        #if os.path.isfile(sOutFile): continue

        sScript  = '%s collapse_elip_by_pos_bychr %s %s %s %s %s %s %s '        \
                    % (sTempScript, sChrID, sInDir, sChrDir, sSizeDir, sOutFile, bStrictPeaks, bAllRBPs)

        if bTestRun:
            print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s'
                      % (sScript, sLogDir, sQueue, sJobName, sChrID))
        #if END:
    #loop END: sChrID
#def END: qsub_collase_eclip_data


def collapse_elip_by_pos_bychr (sChrID, sInDir, sChrDir, sSizeDir, sOutFile, bStrictPeaks, bAllRBPs):
    print(green('Collapse Eclip by Position and Chr', 1))

    bStrictPeaks   = int(bStrictPeaks)
    bAllRBPs       = int(bAllRBPs)
    list_cClipData = get_clipdata_list(bAllRBPs)

    dict_sRBPID    = {}
    for cClipData in list_cClipData:
        if cClipData.sGeneSym not in dict_sRBPID:
            dict_sRBPID[cClipData.sGeneSym] = []
        dict_sRBPID[cClipData.sGeneSym].append(cClipData.sAccID)
    #loop END: cClipData

    dict_nPos      = {sChrID: {}}
    dict_sOutbyRBP = {sRBPID : [] for sRBPID in dict_sRBPID}
    for sRBPID in dict_sRBPID:

        for sAccID in dict_sRBPID[sRBPID]:
            InFile                  = open('%s/%s.%s.peakdata' % (sInDir, sChrID, sAccID), 'rb')
            list_cPeakData          = pickle.load(InFile)
            InFile.close ()

            dict_sOutbyRBP[sRBPID] += list_cPeakData

            for cPeak in list_cPeakData:

                # Peak Calling Score, 1000 = strict, 200 = lenient
                if bStrictPeaks:
                    if cPeak.nScore != 1000: continue

                sKey = '%s,%s' % (cPeak.nStartPos, cPeak.nEndPos)

                if sKey not in dict_nPos[cPeak.sChrID]:
                    dict_nPos[cPeak.sChrID][sKey] = []

                dict_nPos[cPeak.sChrID][sKey].append([cPeak.sGeneSym, cPeak.sLocale])
            #loop END: cPeak
        #loop END: sAccID
    #loop END: sRBPID

    OutFile = open(sOutFile, 'wb')
    pickle.dump(dict_nPos, OutFile)
    OutFile.close()

    #OutbyRBP
    for sRBPID in dict_sOutbyRBP:

        dict_nSizeRBP = {'intergenic': 0,
                         'intronic':   0,
                         '5UTR':       0,
                         '3UTR':       0,
                         'exonic':     0}

        sOutDir     = '%s/%s' % (sChrDir,  sRBPID)
        sSizeOut    = '%s/%s' % (sSizeDir, sRBPID)
        os.makedirs(sOutDir, exist_ok=True)
        os.makedirs(sSizeOut, exist_ok=True)

        sOutFile    = '%s/peaks_%s_%s.dat'    % (sOutDir, sSEQDATA, sChrID)
        sSizeFile   = '%s/peaksize_%s_%s.dat' % (sSizeOut, sSEQDATA, sChrID)
        list_cPeaks = dict_sOutbyRBP[sRBPID]

        print(sRBPID, len(list_cPeaks))

        dict_nPos = {sChrID: {}}

        for cPeak in list_cPeaks:

            sKey = '%s,%s' % (cPeak.nStartPos, cPeak.nEndPos)

            if sKey not in dict_nPos[cPeak.sChrID]:
                dict_nPos[cPeak.sChrID][sKey] = []

            dict_nPos[cPeak.sChrID][sKey].append([cPeak.sGeneSym, cPeak.sLocale])
            try: dict_nSizeRBP[cPeak.sLocale] += (cPeak.nEndPos - cPeak.nStartPos)
            except KeyError: continue
        #loop END: cPeak

        OutFile = open(sOutFile, 'wb')
        pickle.dump(dict_nPos, OutFile)
        OutFile.close()

        OutFile = open(sSizeFile, 'wb')
        pickle.dump(dict_nSizeRBP, OutFile)
        OutFile.close()
    #loop END: sRBPID
#def END: collapse_elip_by_pos_bychr


def load_eclip_data (sWorkDir, sPrevDir):
    print(green('Load eCLIP data', 1))

    sCLIPFile    = '%s/peak_%s.dat'              % (sWorkDir, sSEQDATA)
    sJustPosFile = '%s/peak_%s_justpos.dat'      % (sWorkDir, sSEQDATA)
    sIDPosFile   = '%s/peak_%s_byid.justpos.dat' % (sWorkDir, sSEQDATA)

    if not os.path.isfile(sJustPosFile):
        if not os.path.isfile(sCLIPFile):
            dict_sOutput = {}

            '''Load CLIP data'''
            list_cClipData = get_clipdata_list()
            sClipDir       = '%s/clipdata/peaks' % sPrevDir

            for cClipData in list_cClipData:

                print('Loaded', cClipData.sAccID)

                InFile         = open('%s/%s.peakdata' % (sClipDir, cClipData.sAccID), 'rb')
                list_cPeakData = pickle.load(InFile)
                InFile.close ()

                if cClipData.sGeneSym not in dict_sOutput:
                    dict_sOutput[cClipData.sGeneSym] = []
                dict_sOutput[cClipData.sGeneSym] += list_cPeakData
            #loop END: cClipData

            #VS-Check
            if not dict_sOutput:
                sys.exit('Invalid List : load_eclip_data : dict_sOutput size= %d' % len(dict_sOutput))

            OutFile = open(sCLIPFile, 'wb')
            pickle.dump(dict_sOutput, OutFile)
            OutFile.close()

            for sGeneSym in dict_sOutput:
                sIndyRBPFile = '%s/%s.dat' % (sClipDir, sGeneSym)

                OutFile = open(sIndyRBPFile, 'wb')
                pickle.dump(dict_sOutput[sGeneSym], OutFile)
                OutFile.close()
            #loop END: sGeneSym

        else:
            InFile       = open(sCLIPFile, 'rb')
            dict_sOutput = pickle.load(InFile)
            InFile.close()
        #if END:

        dict_nPos   = {}
        dict_sID    = {}
        for sAccID in dict_sOutput:

            print('Processed', sAccID)

            list_cPeak = dict_sOutput[sAccID]

            if sAccID not in dict_sID:
                dict_sID[sAccID] = {}


            for cPeak in list_cPeak:

                sKey = '%s,%s' % (cPeak.nStartPos, cPeak.nEndPos)

                if cPeak.sChrID not in dict_nPos:
                    dict_nPos[cPeak.sChrID] = {}

                if cPeak.sChrID not in dict_sID[sAccID]:
                    dict_sID[sAccID][cPeak.sChrID] = {}

                if sKey not in dict_nPos[cPeak.sChrID]:
                    dict_nPos[cPeak.sChrID][sKey] = []

                if sKey not in dict_sID[sAccID][cPeak.sChrID]:
                    dict_sID[sAccID][cPeak.sChrID][sKey] = []

                dict_nPos[cPeak.sChrID][sKey].append([sAccID, cPeak.sLocale])
                dict_sID[sAccID][cPeak.sChrID][sKey].append([cPeak.sGeneSym, cPeak.sLocale])
            #loop END: cPeak
        #loop END: sRBP

        OutFile = open(sJustPosFile, 'wb')
        pickle.dump(dict_nPos, OutFile)
        OutFile.close()

        OutFile = open(sIDPosFile, 'wb')
        pickle.dump(dict_sID, OutFile)
        OutFile.close()
    else:
        InFile      = open(sJustPosFile, 'rb')
        dict_nPos   = pickle.load(InFile)
        InFile.close()

        InFile      = open(sIDPosFile, 'rb')
        dict_sID    = pickle.load(InFile)
        InFile.close()
    #if END:
    return dict_nPos, dict_sID
#def END: load_eclip_data


def load_mirna_data (sPrevDir):
    print(green('Load miR data', 1))
    #Load miR target file
    sInFile         = '%s/mir_targets_top%s%%_cap.dat' % (sPrevDir, int(fTOP_EXPR*100))
    InFile          = open(sInFile, 'rb')
    dict_sTargetPos = pickle.load(InFile)
    InFile.close()

    list_nSize      = []
    for sChrID in dict_sTargetPos:

        for sGeneSym in dict_sTargetPos[sChrID]:

            for sPosKey in dict_sTargetPos[sChrID][sGeneSym]:

                nStartPos, nEndPos = [int(sPos) for sPos in sPosKey.split('-')]
                list_nSize.append(nEndPos - nStartPos)
            #loop ENd sPosKey
        #loop ENd sGeneSym
    #loop ENd sChrID

    print('MirTar Size', sum(list_nSize))
    return dict_sTargetPos
#def END: load_mirna_data


def load_refseq_data (sPrevDir):
    print(green('Load RefSeq data', 1))
    if sSEQDATA == 'WXS':
        sCaptureName  = sCAPTUREFILE.split('/')[-1].split('.')[0]
        sInFile       = '%s/ref_capture/refseq_%s.dat' % (sPrevDir, sCaptureName)
    else:
        sInFile       = '%s/refseq.dat' % sPrevDir

    print(sInFile)

    InFile       = open(sInFile, 'rb')
    list_cRefSeq = pickle.load(InFile)
    InFile.close()

    if sSEQDATA == 'WGS':
        return [list_cRefSeq[sKey] for sKey in list_cRefSeq]
    else:
        return list_cRefSeq
#def load_refseq_data


def qsub_assign_variants (sWorkDir, sClipDir, sMirDir, nBufferSize, nChrWin, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun):

    sRBPGroup   = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall   = 'strict' if bStrictPeaks else 'allpeaks'
    sVCFDir     = '%s/VCFs-%s/VCFs_bychrwin'            % (sWorkDir, sSEQDATA)
    sChrDir     = '%s/%s/%s-RBPs_bychr'                 % (sClipDir, sPeakCall, sRBPGroup)
    sOutDir     = '%s/variant_assigned_bychr_%s_%s_%s'  % (sWorkDir, sPeakCall, sRBPGroup, sSEQDATA)

    sPrevJob    = 'Jinman.CollapseByPos.VCF.%s'         % sSTUDY
    sJobName    = 'Jinman.VarAssignmentByChr.%s.%s.%s'  % (sSEQDATA, sPeakCall, sRBPGroup)
    sLogDir     = '%s/log/%s/%s/%s'                     % (sBASE_DIR, sJobName, sSTUDY, sTIME_STAMP)
    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    for sChrID in list_sCHRIDs:

        list_nWindow    = set_chr_range(sChrID, nChrWin)

        for i, sViewRange in enumerate(list_nWindow):

            sInterval = '%s.%s'                 % (i, sViewRange)
            sVarFile  = '%s/vcfs_%s_%s.%s.dat'  % (sVCFDir, sSEQDATA, sChrID, sViewRange)
            sPeakFile = '%s/peaks_%s_%s.dat'    % (sChrDir, sSEQDATA, sChrID)

            sOutFile  = '%s/assigned_variants_%s.%s.dat' % (sOutDir, sChrID, sViewRange)

            sScript   = '%s assign_variants_bychr %s %s %s %s %s %s'        \
                        % (sTempScript, sChrID, sVarFile, sPeakFile, sMirDir, sOutFile, nBufferSize)

            if bTestRun:
                print(sScript)
            else:
                os.makedirs(sLogDir, exist_ok=True)
                os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s.%s -hold_jid %s.%s.*'
                          % (sScript, sLogDir, sQueue, sJobName, sChrID, sInterval, sPrevJob, sChrID))
        #loop END: i, sViewRange
    #loop END: sChrID
#def END: qsub_assign_variants


def qsub_assign_variants_byRBP (sWorkDir, sClipDir, sMirDir, nBufferSize, nChrWin, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun):

    sRBPGroup   = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall   = 'strict' if bStrictPeaks else 'allpeaks'
    sVCFDir     = '%s/VCFs-%s/VCFs_bychrwin'                 % (sWorkDir, sSEQDATA)
    sChrDir     = '%s/%s/%s-RBPs_bychr'                      % (sClipDir, sPeakCall, sRBPGroup)
    sOutDir     = '%s/variant_assigned_byRBP_%s_%s_%s'       % (sWorkDir, sPeakCall, sRBPGroup, sSEQDATA)
    sJobName    = 'Jinman.VarAssignmentByRBPChr.%s.%s.%s.%s' % (sSTUDY, sPeakCall, sRBPGroup, sSEQDATA)
    sLogDir     = '%s/log/%s/%s'                             % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    for sChrID in list_sCHRIDs:

        if sChrID != 'chr11': continue

        list_nWindow    = set_chr_range(sChrID, nChrWin)

        for i, sViewRange in enumerate(list_nWindow):

            sInterval = '%s.%s'                          % (i, sViewRange)
            sVarFile  = '%s/vcfs_%s_%s.%s.dat'           % (sVCFDir, sSEQDATA, sChrID, sViewRange)

            sScript   = '%s assign_variants_byrbp %s %s %s %s %s %s %s %s'        \
                        % (sTempScript, sChrID, sVarFile, sOutDir, sChrDir, sMirDir, sViewRange, bAllRBPs, nBufferSize)

            if bTestRun:
                print(sScript)
            else:
                os.makedirs(sLogDir, exist_ok=True)
                os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s.%s'
                          % (sScript, sLogDir, sQueue, sJobName, sChrID, sInterval))
        #loop END: i, sViewRange
    #loop END: sChrID
#def END: qsub_assign_variants_byRBP


def qsub_combine_variants_byRBP (sWorkDir, nChrWin, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun):

    sRBPGroup   = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall   = 'strict' if bStrictPeaks else 'allpeaks'
    sVarDir     = '%s/variant_assigned_byRBP_%s_%s_%s'             % (sWorkDir, sPeakCall, sRBPGroup, sSEQDATA)
    sPrevJob    = 'Jinman.VarAssignmentByRBPChr.%s.%s.%s.%s'       % (sSTUDY, sPeakCall, sRBPGroup, sSEQDATA)
    sJobName    = 'Jinman.Combine.VarAssignmentByRBPChr.%s.%s.%s'  % (sSTUDY, sPeakCall, sRBPGroup)
    sLogDir     = '%s/log/%s/%s'                                   % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sLogDir, exist_ok=True)

    list_cClipData = get_clipdata_list(bAllRBPs)

    dict_sRBPID    = {}
    for cClipData in list_cClipData:
        if cClipData.sGeneSym not in dict_sRBPID:
            dict_sRBPID[cClipData.sGeneSym] = []
        dict_sRBPID[cClipData.sGeneSym].append(cClipData.sAccID)
    #loop END: cClipData

    list_sRBPID = list(dict_sRBPID.keys())

    for sRBPID in list_sRBPID:

        #sOutFile  = '%s/%s_variant_assigned_byRBP_%s.dat' % (sVarDir, sRBPID, sSEQDATA)
        #sScript   = '%s combine_variant_assignment %s %s %s %s' \
        #             % (sTempScript, sVarDir, sRBPID, nChrWin, sOutFile)

        # Just RBPOnly *Examine Top RBP ###
        sOutFile = '%s/%s_variant_assigned_byRBP_%s-RBPOnly.dat' % (sVarDir, sRBPID, sSEQDATA)
        sScript   = '%s combine_variant_assignment_justRBP %s %s %s %s' \
                    % (sTempScript, sVarDir, sRBPID, nChrWin, sOutFile)
        ###################################
        if bTestRun:
            print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s -hold_jid %s.*'
                      % (sScript, sLogDir, sQueue, sJobName, sRBPID, sPrevJob))
    #loop END: sRBPID
#def END: qsub_combine_variants_byRBP


def assign_variants_bychr (sChrID, sVarFile, sPeakFile, sMirDir, sOutFile, nBufferSize):
    print(green('Assign variants by chr', 1))
    list_sRegion   = ['intergenic', 'intronic', 'exonic', '5UTR', '3UTR', 'other']
    list_sRBP      = []
    list_sMir      = []
    list_nAllPos   = []
    list_nForInter = []
    nBufferSize    = int(nBufferSize)

    #Load miRNA data
    dict_sMirs   = load_mirna_data (sMirDir)
    print('dict_sMirs', len(dict_sMirs))

    InFile         = open(sPeakFile, 'rb')
    dict_sPeakPos  = pickle.load(InFile)
    InFile.close()
    print('dict_sPeakPos', len(dict_sPeakPos))

    InFile         = open(sVarFile, 'rb')
    dict_sVarPos   = pickle.load(InFile)
    InFile.close()
    print('dict_sVarPos', len(dict_sVarPos))

    dict_nPeakPos = bin_peak_positions (dict_sPeakPos[sChrID], nBufferSize, 0)

    for nPos in dict_sVarPos:
        bMatch  = 0
        sKey    = nPos

        list_nAllPos.append(sKey)

        if check_rbp (nPos, dict_nPeakPos):
            bMatch = 1
            list_sRBP.append(sKey)

        if check_mir (nPos, dict_sMirs[sChrID], nBufferSize, 0):
            bMatch = 1
            list_sMir.append(sKey)

        if not bMatch:
            list_nForInter.append(sKey)
    #loop END: nPos

    # Both Mir and RBP
    list_sRBPMir  = list(set(list_sRBP) & set(list_sMir))

    # RBP Only
    list_sRBPOnly = [nPos for nPos in list_sRBP if nPos not in set(list_sRBPMir)]
    list_sMirOnly = [nPos for nPos in list_sMir if nPos not in set(list_sRBPMir)]

    #list_nInter = list(set(list_nAllPos) - set(list_sRBP + list_sMir))

    list_nInterMatch = []
    list_nNoSite     = []

    dict_nPeakPos    = bin_peak_positions (dict_sPeakPos[sChrID], nBufferSize, 1)
    for nPos in list_nForInter:

        bMatch        = 0
        sKey          = nPos

        list_nAllPos.append(sKey)

        if check_rbp (nPos, dict_nPeakPos):
            bMatch = 1

        if check_mir (nPos, dict_sMirs[sChrID], nBufferSize, 1):
            bMatch = 1

        if bMatch:
            list_nInterMatch.append(sKey)

        if not bMatch:
            list_nNoSite.append(sKey)
    #loop END: nPos

    dict_sClass  = {'RBPOnly':list_sRBPOnly,
                    'MirOnly':list_sMirOnly,
                    'RBPMir' :list_sRBPMir,
                    'Inter'  :list_nInterMatch,
                    'NoSite' :list_nNoSite
                    }

    dict_sOutput = {}
    for sClass in dict_sClass:
        if sClass not in dict_sOutput:
            dict_sOutput[sClass] = ''

        list_nPos    = dict_sClass[sClass]
        dict_sRegion = {sRegion : {} for sRegion in list_sRegion}

        for nPos in list_nPos:
            sLocale  = list(set([cVCF.sLocale for cVCF in dict_sVarPos[nPos]]))[0]

            if nPos not in dict_sRegion[sLocale]:
                dict_sRegion[sLocale][nPos] = ''
            dict_sRegion[sLocale][nPos] = len(dict_sVarPos[nPos])
        #loop END: nPos
        dict_sOutput[sClass] = dict_sRegion
    #loop END: sClass

    print('All Variants', len(list_nAllPos))
    for sClass in dict_sOutput:
        print(sClass, '\t'.join(['%s' % len(dict_sOutput[sClass][sRegion]) for sRegion in list_sRegion]))

    OutFile  = open(sOutFile, 'wb')
    pickle.dump(dict_sOutput, OutFile)
    OutFile.close()
#def END: assign_variants_bychr


def assign_variants_byrbp(sChrID, sVarFile, sOutDir, sPeakDir, sMirDir, sViewRange, bAllRBPs, nBufferSize):
    print(green('Assign variants by chr', 1))
    list_sRegion   = ['intergenic', 'intronic', 'exonic', '5UTR', '3UTR', 'other']

    nBufferSize    = int(nBufferSize)
    bAllRBPs       = int(bAllRBPs)

    #Load miRNA data
    dict_sMirs     = load_mirna_data (sMirDir)
    print('dict_sMirs', len(dict_sMirs))

    InFile         = open(sVarFile, 'rb')
    dict_sVarPos   = pickle.load(InFile)
    InFile.close()
    print('dict_sVarPos', len(dict_sVarPos))

    list_cClipData = get_clipdata_list(bAllRBPs)

    dict_sRBPID    = {}
    for cClipData in list_cClipData:
        if cClipData.sGeneSym not in dict_sRBPID:
            dict_sRBPID[cClipData.sGeneSym] = []
        dict_sRBPID[cClipData.sGeneSym].append(cClipData.sAccID)
    #loop END: cClipData

    list_sRBPID = list(dict_sRBPID.keys())
    sTargetRBP  = 'HNRNPU-human'

    for sRBPID in list_sRBPID:

        #################################
        if sRBPID != sTargetRBP: continue
        #################################
        list_sRBP       = []
        list_sMir       = []
        list_nAllPos    = []
        list_nForInter  = []
        print(sRBPID)

        sPeakFile       = '%s/%s/peaks_%s_%s.dat'    % (sPeakDir, sRBPID, sSEQDATA, sChrID)
        InFile          = open(sPeakFile, 'rb')
        dict_sPeakPos   = pickle.load(InFile)
        InFile.close()

        sOutDir2        = '%s/%s' % (sOutDir, sRBPID)
        os.makedirs(sOutDir2, exist_ok=True)
        sOutFile        = '%s/assigned_variants_%s.%s.dat' % (sOutDir2, sChrID, sViewRange)

        dict_nPeakPos   = bin_peak_positions (dict_sPeakPos[sChrID], nBufferSize, 0)

        for nPos in dict_sVarPos:
            bMatch  = 0
            sKey    = nPos

            list_nAllPos.append(sKey)

            if check_rbp (nPos, dict_nPeakPos):
                bMatch = 1
                list_sRBP.append(sKey)

            if check_mir (nPos, dict_sMirs[sChrID], nBufferSize, 0):
                bMatch = 1
                list_sMir.append(sKey)

            if not bMatch:
                list_nForInter.append(sKey)
        #loop END: nPos

        # Both Mir and RBP
        list_sRBPMir    = list(set(list_sRBP) & set(list_sMir))

        # RBP Only
        list_sRBPOnly   = [nPos for nPos in list_sRBP if nPos not in set(list_sRBPMir)]

        # Mir Only
        list_sMirOnly   = [nPos for nPos in list_sMir if nPos not in set(list_sRBPMir)]


        #list_nInter = list(set(list_nAllPos) - set(list_sRBP + list_sMir))

        list_nInterMatch = []
        list_nNoSite     = []

        dict_nPeakPos    = bin_peak_positions (dict_sPeakPos[sChrID], nBufferSize, 1)
        for nPos in list_nForInter:

            bMatch        = 0
            sKey          = nPos

            list_nAllPos.append(sKey)

            if check_rbp (nPos, dict_nPeakPos):
                bMatch = 1

            if check_mir (nPos, dict_sMirs[sChrID], nBufferSize, 1):
                bMatch = 1

            if bMatch:
                list_nInterMatch.append(sKey)

            if not bMatch:
                list_nNoSite.append(sKey)
        #loop END: nPos

        dict_sClass  = {'RBPOnly':list_sRBPOnly ,
                        'MirOnly':list_sMirOnly ,
                        'RBPMir':list_sRBPMir ,
                        'Inter':list_nInterMatch ,
                        'NoSite':list_nNoSite}

        dict_sOutput = {}
        for sClass in dict_sClass:

            if sClass != 'RBPOnly': continue

            if sClass not in dict_sOutput:
                dict_sOutput[sClass] = ''

            list_nPos    = dict_sClass[sClass]
            dict_sRegion = {sRegion : {} for sRegion in list_sRegion}

            for nPos in list_nPos:
                print(nPos, [cVCF.sLocale for cVCF in dict_sVarPos[nPos]])

                sLocale = list(set([cVCF.sLocale for cVCF in dict_sVarPos[nPos]]))[0]

                if nPos not in dict_sRegion[sLocale]:
                    dict_sRegion[sLocale][nPos] = ''
                dict_sRegion[sLocale][nPos] = len(dict_sVarPos[nPos])


            #loop END: nPos
            print(dict_sRegion)

            sys.exit()

            dict_sOutput[sClass] = dict_sRegion
        #loop END: sClass

        print('All Variants', len(list_nAllPos))
        for sClass in dict_sOutput:
            print(sClass, '\t'.join(['%s' % len(dict_sOutput[sClass][sRegion]) for sRegion in list_sRegion]))

        sys.exit()


        OutFile  = open(sOutFile, 'wb')
        pickle.dump(dict_sOutput, OutFile)
        OutFile.close()
    #loop END: sRBPID
#def END: assign_variants_byrbp


def combine_variant_assignment (sVarDir, sRBPID, nChrWin, sOutFile):

    list_sOutputKey = ['RBPOnly', 'MirOnly', 'RBPMir', 'Inter', 'NoSite']
    list_sRegion    = ['intergenic', 'intronic', 'exonic', '5UTR', '3UTR', 'other']

    nChrWin         = int(nChrWin)
    sInDir          = '%s/%s' % (sVarDir, sRBPID)

    dict_sOutput    = {sRegion: {sClass:{} for sClass in list_sOutputKey} for sRegion in list_sRegion}
    for sChrID in list_sCHRIDs:

        list_nWindow  = set_chr_range(sChrID, nChrWin)

        for i, sViewRange in enumerate(list_nWindow):

            sInFile   = '%s/assigned_variants_%s.%s.dat' % (sInDir, sChrID, sViewRange)
            InFile    = open(sInFile, 'rb')
            dict_nVar = pickle.load(InFile)
            InFile.close()

            for sClass in list_sOutputKey:
                for sRegion in list_sRegion:
                    for nPos in dict_nVar[sClass][sRegion]:
                        nPatCnt = dict_nVar[sClass][sRegion][nPos]
                        if nPos not in dict_sOutput[sRegion][sClass]:
                            dict_sOutput[sRegion][sClass][nPos] = 0
                        dict_sOutput[sRegion][sClass][nPos] = nPatCnt
                    #loop END: nPos
                #loop END: sRegion
            #loop END: sClass
        #loop END: i, sViewRange
    #loop END: sChrID

    OutFile = open(sOutFile, 'wb')
    pickle.dump(dict_sOutput, OutFile)
    OutFile.close()
#def END: combine_variant_assignment


def combine_variant_assignment_justRBP (sVarDir, sRBPID, nChrWin, sOutFile):

    list_sRegion      = ['intergenic', 'intronic', 'exonic', '5UTR', '3UTR', 'other']

    nChrWin           = int(nChrWin)
    sInDir            = '%s/%s' % (sVarDir, sRBPID)

    dict_sOutput      = {sRegion: {} for sRegion in list_sRegion}
    for sChrID in list_sCHRIDs:
        print('Loading', sChrID)
        list_nWindow  = set_chr_range(sChrID, nChrWin)

        for i, sViewRange in enumerate(list_nWindow):

            sInFile   = '%s/assigned_variants_%s.%s.dat' % (sInDir, sChrID, sViewRange)
            InFile    = open(sInFile, 'rb')
            dict_nVar = pickle.load(InFile)
            InFile.close()

            sClass    = 'RBPOnly'

            for sRegion in list_sRegion:
                if sChrID not in dict_sOutput[sRegion]:
                    dict_sOutput[sRegion][sChrID] = {}

                for nPos in dict_nVar[sClass][sRegion]:
                    nPatCnt = dict_nVar[sClass][sRegion][nPos]

                    if nPos not in dict_sOutput[sRegion][sChrID]:
                        dict_sOutput[sRegion][sChrID][nPos] = ''
                    dict_sOutput[sRegion][sChrID][nPos] = nPatCnt
                #loop END: nPos
            #loop END: sRegion
        #loop END: i, sViewRange
    #loop END: sChrID
    '''
    for sRegion in dict_sOutput:
        print(sRegion, len(dict_sOutput[sRegion]))
        for sChrID in dict_sOutput[sRegion]:
            print(sChrID, dict_sOutput[sRegion][sChrID])
    '''
    OutFile = open(sOutFile, 'wb')
    pickle.dump(dict_sOutput, OutFile)
    OutFile.close()
#def END: combine_variant_assignment_justRBP



def assign_variants (sWorkDir, dict_cVCF, dict_sMirs, dict_cPeak, dict_sGeneSym, nBufferSize, nVarWindow):
    print(green('Assign variants to location', 1))

    sAssigedFile = '%s/assigned_variants_%s_buffer%s_win%s_top%s.dat' % (sWorkDir, sSEQDATA, nBufferSize, nVarWindow, int(fTOP_EXPR*100))

    if not os.path.isfile(sAssigedFile):

        list_sRBP      = []
        list_sMir      = []
        list_nAllPos   = []
        list_nForInter = []

        print(time.ctime())
        for sChrID in dict_cVCF:

            dict_nPeakPos = bin_peak_positions (dict_cPeak[sChrID], nBufferSize, 0)

            for nPos in dict_cVCF[sChrID]:

                bMatch = 0
                sKey   = '%s,%s' % (sChrID, nPos)
                cVCF   = dict_cVCF[sChrID][nPos][0]
                try: fFPKM  = dict_fFPKM[cVCF.sGeneSym]
                except KeyError: continue

                list_nAllPos.append(sKey)

                if check_rbp (nPos, dict_nPeakPos):
                    bMatch = 1
                    list_sRBP.append(sKey)

                if check_mir (nPos, dict_sMirs[sChrID], nBufferSize, 0):
                    bMatch = 1
                    list_sMir.append(sKey)

                if not bMatch:
                    list_nForInter.append(sKey)

            #loop END: nPos
        #loop END: sChrID

        # Both Mir and RBP
        list_sRBPMir  = list(set(list_sRBP) & set(list_sMir))

        # RBP Only
        list_sRBPOnly = [nPos for nPos in list_sRBP if nPos not in set(list_sRBPMir)]
        list_sMirOnly = [nPos for nPos in list_sMir if nPos not in set(list_sRBPMir)]

        print('All Variants', len(list_nAllPos))
        print('Found RBP Match', len(list_sRBP))
        print('Found Mir Match', len(list_sMir))
        print('RBP-Mir Match', len(list_sRBPMir))
        print('RBP only', len(list_sRBPOnly))
        print('Mir only', len(list_sMirOnly))

        #list_nInter = list(set(list_nAllPos) - set(list_sRBP + list_sMir))

        print('For Intermediate Survey', len(list_nForInter))

        list_nInterMatch = []
        list_nNoSite     = []
        dict_nForInter   = {}
        for sPosKey in list_nForInter:

            sChrID        = sPosKey.split(',')[0]
            nPos          = int(sPosKey.split(',')[1])
            if sChrID not in dict_nForInter:
                dict_nForInter[sChrID] = []
            dict_nForInter[sChrID].append(nPos)
        #loop END: sPosKey

        for sChrID in dict_nForInter:

            #print(sChrID, len(dict_cVCF[sChrID]))

            dict_nPeakPos = bin_peak_positions (dict_cPeak[sChrID], nBufferSize, 1)

            for nPos in dict_cVCF[sChrID]:

                bMatch        = 0
                sKey          = '%s,%s' % (sChrID, nPos)

                list_nAllPos.append(sKey)

                if check_rbp (nPos, dict_nPeakPos):
                    bMatch = 1

                if check_mir (nPos, dict_sMirs[sChrID], nBufferSize, 1):
                    bMatch = 1

                if bMatch:
                    list_nInterMatch.append(sKey)

                if not bMatch:
                    list_nNoSite.append(sKey)
            #loop END: nPos
        #loop END: sPosKey

        print('Intermediate', len(list_nInterMatch))
        print('NoSite', len(list_nNoSite))
        print(time.ctime())

        dict_sOutput = {'RBPOnly':list_sRBPOnly ,
                        'MirOnly':list_sMirOnly ,
                        'RBPMir':list_sRBPMir ,
                        'Inter':list_nInterMatch ,
                        'NoSite':list_nNoSite
                        }

        OutFile  = open(sAssigedFile, 'wb')
        pickle.dump(dict_sOutput, OutFile)
        OutFile.close()
    else:
        InFile       = open(sAssigedFile, 'rb')
        dict_sOutput = pickle.load(InFile)
        InFile.close()

        list_nAllPos = []
        for sClass in dict_sOutput:
            list_nAllPos += dict_sOutput[sClass]
            print(sClass, len(dict_sOutput[sClass]))
        #loop END:
        print('All Pos', len(list_nAllPos))

    return dict_sOutput
#def END: assign_variants


def assign_variants_by_rbp (sWorkDir, dict_cVCF, dict_sMirs, dict_cPeak, dict_fFPKM, nBufferSize, nVarWindow):
    print(green('Assign variants to location', 1))

    sOutDir = '%s/clipdata/byRBPs' % sWorkDir
    os.makedirs(sOutDir, exist_ok=True)

    for sRBP in dict_cPeak:

        sAssigedFile   = '%s/%s_assigned_variants_buffer%s_win%s_top%s.dat' \
                          % (sOutDir, sRBP, nBufferSize, nVarWindow, int(fTOP_EXPR*100))

        if os.path.isfile(sAssigedFile): continue
        print(sRBP, len(dict_cPeak[sRBP]))

        list_sRBP      = []
        list_sMir      = []
        list_nAllPos   = []
        list_nForInter = []

        print(time.ctime())
        for sChrID in dict_cVCF:

            try: dict_cPeak[sRBP][sChrID]
            except KeyError: continue

            dict_nPeakPos = bin_peak_positions (dict_cPeak[sRBP][sChrID], nBufferSize, 0)

            for nPos in dict_cVCF[sChrID]:

                bMatch = 0
                sKey   = '%s,%s' % (sChrID, nPos)
                cVCF   = dict_cVCF[sChrID][nPos][0]
                try: fFPKM  = dict_fFPKM[cVCF.sGeneSym]
                except KeyError: continue

                list_nAllPos.append(sKey)

                if check_rbp (nPos, dict_nPeakPos):
                    bMatch = 1
                    list_sRBP.append(sKey)

                if check_mir (nPos, dict_sMirs[sChrID], nBufferSize, 0):
                    bMatch = 1
                    list_sMir.append(sKey)

                if not bMatch:
                    list_nForInter.append(sKey)

            #loop END: nPos
        #loop END: sChrID

        # Both Mir and RBP
        list_sRBPMir  = list(set(list_sRBP) & set(list_sMir))

        # RBP Only
        list_sRBPOnly = [nPos for nPos in list_sRBP if nPos not in set(list_sRBPMir)]
        list_sMirOnly = [nPos for nPos in list_sMir if nPos not in set(list_sRBPMir)]

        print('All Variants', len(list_nAllPos))
        print('Found RBP Match', len(list_sRBP))
        print('Found Mir Match', len(list_sMir))
        print('RBP-Mir Match', len(list_sRBPMir))
        print('RBP only', len(list_sRBPOnly))
        print('Mir only', len(list_sMirOnly))

        #list_nInter = list(set(list_nAllPos) - set(list_sRBP + list_sMir))

        print('For Intermediate Survey', len(list_nForInter))

        list_nInterMatch = []
        list_nNoSite     = []
        dict_nForInter   = {}
        for sPosKey in list_nForInter:

            sChrID        = sPosKey.split(',')[0]
            nPos          = int(sPosKey.split(',')[1])
            if sChrID not in dict_nForInter:
                dict_nForInter[sChrID] = []
            dict_nForInter[sChrID].append(nPos)
        #loop END: sPosKey

        for sChrID in dict_nForInter:

            #print(sChrID, len(dict_cVCF[sChrID]))

            dict_nPeakPos = bin_peak_positions (dict_cPeak[sRBP][sChrID], nBufferSize, 1)

            for nPos in dict_cVCF[sChrID]:

                bMatch        = 0
                sKey          = '%s,%s' % (sChrID, nPos)

                list_nAllPos.append(sKey)

                if check_rbp (nPos, dict_nPeakPos):
                    bMatch = 1

                if check_mir (nPos, dict_sMirs[sChrID], nBufferSize, 1):
                    bMatch = 1

                if bMatch:
                    list_nInterMatch.append(sKey)

                if not bMatch:
                    list_nNoSite.append(sKey)
            #loop END: nPos
        #loop END: sPosKey

        print('Intermediate', len(list_nInterMatch))
        print('NoSite', len(list_nNoSite))
        print(time.ctime())

        dict_sOutput = {'RBPOnly':list_sRBPOnly ,
                        'MirOnly':list_sMirOnly ,
                        'RBPMir':list_sRBPMir ,
                        'Inter':list_nInterMatch ,
                        'NoSite':list_nNoSite
                        }

        OutFile  = open(sAssigedFile, 'wb')
        pickle.dump(dict_sOutput, OutFile)
        OutFile.close()
    #loop END: sRBP
#def END: assign_variants_by_rbp


def bin_peak_positions (list_sPeakPos, nBufferSize, bBuffered):

    list_nPeakPos = [[int(sPosKey.split(',')[0]), int(sPosKey.split(',')[1])] for sPosKey in list_sPeakPos]
    if bBuffered:
        list_nPeakPos = [[nStart-nBufferSize, nEnd+nBufferSize] for nStart, nEnd in list_nPeakPos]

    list_nPeakPos = sorted(list_nPeakPos, key=lambda e:e[0])

    nBins         = 20
    dict_nBin     = {i:[] for i in range(nBins)}
    nTotalSize    = len(list_nPeakPos)

    for i, (sStart, nEnd) in enumerate(list_nPeakPos):
        nBin = int(i/nTotalSize * nBins)
        dict_nBin[nBin].append([sStart, nEnd])
    #loop END:

    dict_nPos     = {}
    for nBin in dict_nBin:
        if not dict_nBin[nBin]: continue
        nBinStart = min([nStart for nStart,nEnd in dict_nBin[nBin]])
        nBinEnd   = max([nEnd   for nStart,nEnd in dict_nBin[nBin]])

        sKey      = '%s,%s' % (nBinStart, nBinEnd)
        if sKey not in dict_nPos:
            dict_nPos[sKey] = []
        dict_nPos[sKey] = dict_nBin[nBin]
    #loop END: nBin

    return dict_nPos
#def END: bin_peak_positions


def check_rbp (nPos, dict_nPos):

    bFlag         = 0
    for sBinKey in dict_nPos:
        nBinStart, nBinEnd = [int(sPos) for sPos in sBinKey.split(',')]

        if not (nBinStart <= nPos <= nBinEnd): continue

        for nStartPos, nEndPos in dict_nPos[sBinKey]:
            if nStartPos <= nPos <= nEndPos:
                bFlag = 1
                break
        #loop END: nStartPos, nEndPos
    return bFlag
#def END: check_rbp


def check_mir (nPos, dict_sTargets, nBufferSize, bBuffered):
    bFlag = 0
    for sGeneSym in dict_sTargets:

        for sPosKey in dict_sTargets[sGeneSym]:

            nStartPos, nEndPos = [int(sPos) for sPos in sPosKey.split('-')]

            if bBuffered:
                nStartPos = nStartPos - nBufferSize
                nEndPos   = nEndPos   + nBufferSize


            if nStartPos <= nPos <= nEndPos:
                bFlag = 1
                break
        #loop END: sPoskey
    #loop END: sGeneSym
    return bFlag
#def END: check_mir


def calculate_regionsize (sWorkDir, sPrevDir, sPrev2Dir, list_cRefSeq, bStrictPeaks, bAllRBPs):
    print(green('Calculate region sizes', 1))

    sRBPGroup        = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall        = 'strict' if bStrictPeaks else 'allpeaks'

    list_sRegion     = ['intergenic', 'intronic', 'exonic', '5UTR', '3UTR']
    list_sRegion     = ['3UTR', '5UTR', 'exonic', 'intergenic', 'intronic']

    #dict_nLockedSize   = {'intergenic':     4810, 'intronic':     946816,
    #                      '5UTR':         397003, '3UTR':         283951,
    #                      'exonic':      1900873}
    dict_nRefSize_filt = calculate_filtered_refseqlen (list_cRefSeq, list_sRegion)
    dict_nRBPSize_filt = calculate_filtered_rbplen    (sPrevDir, list_sRegion)
    dict_nRefSize_full = calculate_full_refseqlen     (sPrev2Dir, list_sRegion)
    dict_nRBPSize_full = calculate_full_rbplen        (sPrev2Dir, list_sRegion, bStrictPeaks, bAllRBPs)

    for sRegion in list_sRegion:

        sOut = '%s\t%s\t%s\t%s\t%s\n' % (sRegion, dict_nRBPSize_full[sRegion], dict_nRefSize_full[sRegion],
                                             dict_nRBPSize_filt[sRegion], dict_nRefSize_filt[sRegion])
        #print(sOut[:-1])
        #OutFile.write(sOut)

        #if sSTUDY == 'TCGA_LIHC':
            #dict_nRefSize_filt[sRegion] = dict_nRefSize_filt[sRegion] - dict_nLockedSize[sRegion]
        #elif sSTUDY == 'TCGA_LAML':
            #dict_nRefSize_filt[sRegion] = dict_nRefSize_filt[sRegion] - dict_nRBPSize_filt[sRegion]
    #loop END: sRegion
    #OutFile.close()

    sOutFile = '%s/%s_RBP_Region_%s_for%s.data' % (sWorkDir, sRBPGroup, sPeakCall, sSEQDATA)
    OutFile  = open(sOutFile, 'wb')

    if sSEQDATA == 'WXS':
        pickle.dump(dict_nRefSize_filt, OutFile)
    else:
        pickle.dump(dict_nRefSize_full, OutFile)
    #if END:
    OutFile.close()
#def END: calculate_regionsize


def calculate_filtered_refseqlen (list_cRefSeq, list_sRegion):

    dict_sOutput    = {sRegion : 0 for sRegion in list_sRegion}

    if sSEQDATA == 'WGS': return dict_sOutput
    else:

        for cRef in list_cRefSeq:
            dict_sOutput['5UTR']       += len(cRef.s5UTRSeq)
            dict_sOutput['3UTR']       += len(cRef.s3UTRSeq)
            dict_sOutput['exonic']     += len(cRef.sORFSeq)
            dict_sOutput['intergenic'] += sum([nEnd - nStart for nStart, nEnd in cRef.list_nIntergenic])
            dict_sOutput['intronic']   += sum([nEnd - nStart for nStart, nEnd in cRef.list_nIntrons])
        #loop END: cRef

        return dict_sOutput
#def END: calculate_filtered_refseq

def calculate_filtered_rbplen (sPrevDir, list_sRegion):

    dict_sOutput    = {sRegion : 0 for sRegion in list_sRegion}
    if sSEQDATA == 'WGS': return dict_sOutput
    else:
        sFilteredPeaks  = '%s/clipdata/%s_combined.sorted.merged.filter.dat' % (sPrevDir, sCELLS)
        InFile        = open(sFilteredPeaks, 'rb')
        list_cPeak    = pickle.load(InFile)
        InFile.close()
        print('list_cPeak', len(list_cPeak))
        for cPeak in list_cPeak:
            if not cPeak.sLocale: continue
            dict_sOutput[cPeak.sLocale] += cPeak.nPeakSize
        #loop END: cPeak

        return dict_sOutput
#def END: calculate_filtered_rbplen


def calculate_full_refseqlen (sPrev2Dir, list_sRegion):

    dict_sOutput = {sRegion : 0 for sRegion in list_sRegion}
    sInFile      = '%s/refseq.dat' % sPrev2Dir
    InFile       = open(sInFile, 'rb')
    dict_cRefSeq = pickle.load(InFile)
    InFile.close()

    for sGeneSym in dict_cRefSeq:
        cRef            = dict_cRefSeq[sGeneSym]
        dict_sOutput['5UTR']      += len(cRef.s5UTRSeq)
        dict_sOutput['3UTR']      += len(cRef.s3UTRSeq)
        dict_sOutput['exonic']    += len(cRef.sORFSeq)
        dict_sOutput['intronic']  += cRef.nIntronSize
    #loop END: sNMID
    dict_sOutput['intergenic'] = nGENOME_SIZE - sum([dict_sOutput[sRegion] for sRegion in dict_sOutput if sRegion != 'Intergenic'])

    return dict_sOutput
#def END: calculate_full_refseqlen


def calculate_full_rbplen (sPrev2Dir, list_sRegion, bStrictPeaks, bAllRBPs):


    sRBPGroup       = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall       = 'strict' if bStrictPeaks else 'allpeaks'

    dict_sOutput    = {sRegion : 0 for sRegion in list_sRegion}

    sFilteredPeaks  = '%s/clipdata/%s/%s_combined.sorted.merged.dat' % (sPrev2Dir, sPeakCall, sRBPGroup)
    print(sFilteredPeaks)

    InFile        = open(sFilteredPeaks, 'rb')
    list_cPeak    = pickle.load(InFile)
    InFile.close()
    print('list_cPeak', len(list_cPeak))
    for cPeak in list_cPeak:
        try: dict_sOutput[cPeak.sLocale] += cPeak.nPeakSize
        except KeyError: continue
    #loop END: cPeak
    return dict_sOutput
#def END: calculate_full_rbplen


def calculate_regionsize_by_rbp_old (sWorkDir, sPrevDir, dict_fFPKM, nBufferSize, nVarWindow):
    print(green('Calculate region sizes by RBP', 1))

    list_cClipData = get_clipdata_list()
    dict_sRBP      = {cClipData.sGeneSym : [] for cClipData in list_cClipData}

    sOutDir       = '%s/clipdata/byRBPs' % sWorkDir
    os.makedirs(sOutDir, exist_ok=True)

    sMergedPeaks  = '%s/clipdata/%s_combined.sorted.merged.filter.dat' % (sPrevDir, sCELLS)
    InFile        = open(sMergedPeaks, 'rb')
    list_cPeak    = pickle.load(InFile)
    InFile.close()
    print('list_cPeak', len(list_cPeak))

    for cPeak in list_cPeak:
        try: dict_fFPKM[cPeak.sGeneSym]
        except KeyError: continue

        for sSample in cPeak.list_sSample:

            sKey = sSample.split('_')[0]
            sKey = '%s-human' % sKey
            dict_sRBP[sKey].append(cPeak)
    #loop END: cPeak

    list_sRBPs    = list(dict_sRBP.keys())
    for sRBP in list_sRBPs:

        sOutFile = '%s/%s_rbpsize_buffer%s_win%s_top%s.dat' \
                    % (sOutDir, sRBP, nBufferSize, nVarWindow, int(fTOP_EXPR * 100))

        print(sRBP, len(dict_sRBP[sRBP]))

        dict_nSizeRBP = {'intergenic': 0,
                         'intronic':   0,
                         '5UTR':       0,
                         '3UTR':       0,
                         'exonic':     0}
        list_cPeak = dict_sRBP[sRBP]
        for cPeak in list_cPeak:
            dict_nSizeRBP[cPeak.sLocale] += cPeak.nPeakSize
        #loop END: sPosKey

        OutFile  = open(sOutFile, 'wb')
        pickle.dump(dict_nSizeRBP, OutFile)
        OutFile.close()
    #loop ENDL sRBP
#def END: calculate_regionsize_by_rbp


def combined_rbpsize (sWorkDir, bStrictPeaks, bAllRBPs):
    print(green('Calculate region sizes by RBP', 1))

    sRBPGroup      = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall      = 'strict' if bStrictPeaks else 'allpeaks'
    sOutDir        = '%s/%s/%s-RBPsizes'    % (sWorkDir, sPeakCall, sRBPGroup)
    list_cClipData = get_clipdata_list(bAllRBPs)

    dict_sRBPID    = {}
    for cClipData in list_cClipData:
        if cClipData.sGeneSym not in dict_sRBPID:
            dict_sRBPID[cClipData.sGeneSym] = []
        dict_sRBPID[cClipData.sGeneSym].append(cClipData.sAccID)
    #loop END: cClipData

    list_sRBPID = list(dict_sRBPID.keys())
    for sRBP in list_sRBPID:

        sOutFile      = '%s/%s_peaksize_combined.dat' % (sOutDir, sRBP)
        print(sOutFile)
        dict_nSizeRBP = {'intergenic': 0,
                         'intronic':   0,
                         '5UTR':       0,
                         '3UTR':       0,
                         'exonic':     0}

        for sChrID in list_sCHRIDs:
            sPeakFile      = '%s/%s/peaksize_%s_%s.dat' % (sOutDir, sRBP, sSEQDATA, sChrID)
            InFile         = open(sPeakFile, 'rb')
            dict_sPeakSize = pickle.load(InFile)
            InFile.close()

            for sLocale in dict_sPeakSize:
                dict_nSizeRBP[sLocale] += dict_sPeakSize[sLocale]
            #loop END: sPosKey
        #loop END: sChrID

        #print(sRBP, dict_nSizeRBP)

        OutFile  = open(sOutFile, 'wb')
        pickle.dump(dict_nSizeRBP, OutFile)
        OutFile.close()
    #loop ENDL sRBP
#def END: combined_rbpsize


def calculate_genic_size (dict_cRefseq, bCapture):
    dict_sOutput  =  {'intergenic': 0,
                     'intronic':    0,
                     '5UTR':        0,
                     '3UTR':        0,
                     'exonic':      0}

    dict_sChrID   = {}
    for sNMID in dict_cRefseq:
        cRef            = dict_cRefseq[sNMID]
        if cRef.sChrID not in dict_sChrID:
            dict_sChrID[cRef.sChrID] = []
        dict_sChrID[cRef.sChrID].append([cRef.sGeneSym, sNMID])

        #print(cRef.sGeneSym, sNMID, len(cRef.s5UTRSeq), len(cRef.sORFSeq), len(cRef.s3UTRSeq))
        dict_sOutput['5UTR']      += len(cRef.s5UTRSeq)
        dict_sOutput['3UTR']      += len(cRef.s3UTRSeq)
        dict_sOutput['exonic']    += len(cRef.sORFSeq)

        if bCapture:
            dict_sOutput['intergenic'] += sum([nEnd - nStart for nStart, nEnd in cRef.list_nIntergenic])
            dict_sOutput['intronic']  += sum([nEnd - nStart for nStart, nEnd in cRef.list_nIntrons])
        else:
            dict_sOutput['intronic'] += cRef.nIntronSize

    #loop END: sNMID
    if not bCapture:
        dict_sOutput['intergenic'] = nGENOME_SIZE - sum([dict_sOutput[sRegion] for sRegion in dict_sOutput if sRegion != 'Intergenic'])

    #for e in dict_sOutput:
        #print(e, dict_sOutput[e])

    return dict_sOutput

    '''
    list_nIntergenic = []
    for sChrID in dict_sChrID:
        dict_sBedPos  = dict_sCapture[sChrID]

        list_sNonGene = []
        for sNMID, sGeneSym in dict_sChrID[sChrID]:

            list_sNonGene += [sPos for sPos in dict_sBedPos if sNMID not in set(dict_sBedPos[sPos])
                              or sGeneSym not in set(dict_sBedPos[sPos])]
        #loop END:
        list_sNonGene = list(set(list_sNonGene))
        print(sChrID, len(list_sNonGene))

        for sPos in list_sNonGene:
            print(dict_sBedPos[sPos])

        sys.exit()
    '''
#def END: calculate_genic_size


def acquire_top_expressed_genes (list_cRefSeq):
    print(green('Get top expressed genes', 1))

    dict_cRefSeq    = {cRef.sGeneSym:cRef for cRef in list_cRefSeq}

    if fTOP_EXPR == 1: return dict_cRefSeq
    else:
        # Load TCGA FPKM File List
        sFPKMDir        = '%s/%s/FPKMfiles/UQ_norm_recode'    % (sBASE_DIR, sSTUDY)
        sFPKMFile       = '%s/%s/%s_rnaseq-UQ.txt'            % (sBASE_DIR, sSTUDY, sSTUDY)
        dict_sFPKMFiles = load_fpkm_list (sFPKMFile)
        dict_sExpr      = get_expr_data (sFPKMDir, dict_sFPKMFiles)

        #Get Top Expressed
        list_sTopGenes = [[sGeneSym, dict_sExpr[sGeneSym]] for sGeneSym in dict_sExpr if dict_sExpr[sGeneSym] > 0]
        list_sInRefSeq = []
        for sGeneSym, fFPKM in list_sTopGenes:
            try: dict_cRefSeq[sGeneSym]
            except KeyError: continue
            list_sInRefSeq.append([sGeneSym, fFPKM])
        #loop END: sGeneSym, fFPKM
        list_sTopGenes = sorted(list_sInRefSeq, key=lambda e: e[1], reverse=True)[:int(len(list_sInRefSeq) * fTOP_EXPR)]

        dict_sTopGenes = dict(list_sTopGenes)

        return dict_sTopGenes
    #if END:
#def END: acquire_top_expressed_genes


def load_variant_assignments (sWorkDir, nChrWin, list_sOutputKey, list_sRegion, bRecurrency, fRecurrency, bPatCnt, bStrictPeaks, bAllRBPs):

    sRBPGroup    = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall    = 'strict' if bStrictPeaks else 'allpeaks'
    sInDir       = '%s/variant_assigned_bychr_%s_%s_%s' % (sWorkDir, sPeakCall, sRBPGroup, sSEQDATA)
    sAnalysisTag =  '%s%s%s%s'                          % (sSEQDATA, '.recurrent' if bRecurrency else '.nonrecurrent',
                                                           '.posxpat' if bPatCnt else '.bypos', '.%s' % sPeakCall)
    sOutFile     = '%s/%s_variant_assigned_%s.data'     % (sWorkDir, sRBPGroup, sAnalysisTag)

    print(sOutFile)
    dict_sOutput = {sRegion: {sClass:{} for sClass in list_sOutputKey} for sRegion in list_sRegion}

    if not os.path.isfile(sOutFile):
        for sChrID in list_sCHRIDs:

            list_nWindow    = set_chr_range(sChrID, nChrWin)

            for i, sViewRange in enumerate(list_nWindow):

                sInFile   = '%s/assigned_variants_%s.%s.dat' % (sInDir, sChrID, sViewRange)
                InFile    = open(sInFile, 'rb')
                dict_nVar = pickle.load(InFile)
                InFile.close()

                for sClass in list_sOutputKey:
                    for sRegion in list_sRegion:
                        for nPos in dict_nVar[sClass][sRegion]:
                            nPatCnt = dict_nVar[sClass][sRegion][nPos]

                            if bRecurrency:
                                if nPatCnt <= int(nSAMPLE_CNT * fRecurrency): continue

                            if nPos not in dict_sOutput[sRegion][sClass]:
                                dict_sOutput[sRegion][sClass][nPos] = ''

                            dict_sOutput[sRegion][sClass][nPos] = nPatCnt
                        #loop END: nPos
                    #loop END: sRegion
                #loop END: sClass
            #loop END: i, sViewRange
        #loop END: sChrID

        OutFile = open(sOutFile, 'wb')
        pickle.dump(dict_sOutput, OutFile)
        OutFile.close()

    else:
        InFile = open(sOutFile, 'rb')
        dict_sOutput = pickle.load(InFile)
        InFile.close()

    #if END:
    return dict_sOutput
#def END: load_variant_assignments


def normalize_enrich_bychr (sWorkDir, nChrWin, fRecurrency, bPatCnt, bStrictPeaks, bAllRBPs):
    print(green('Determine enrichment by chr %s' % ('PosxPat' if bPatCnt else 'ByPos'), 1))

    bRecurrency      = 0
    list_sOutputKey  = ['RBPOnly', 'MirOnly', 'RBPMir', 'NoSite']
    list_sRegion     = ['intergenic', 'intronic', 'exonic', '5UTR', '3UTR']

    sHeader          = '\t'.join([sClass for sClass in list_sOutputKey])

    print('region\t%s' % sHeader)
    dict_sVars       = load_variant_assignments (sWorkDir, nChrWin, list_sOutputKey, list_sRegion,
                                                 bRecurrency, fRecurrency, bPatCnt, bStrictPeaks, bAllRBPs)

    for sRegion in list_sRegion:
        if bPatCnt:
            list_nPatCnt = [sum([dict_sVars[sRegion][sClass][nPos]
                                 for nPos in dict_sVars[sRegion][sClass]]) for sClass in list_sOutputKey]
            sOut    = '\t'.join(['%s' % nPatCnt for nPatCnt in list_nPatCnt])
        else:
            sOut    = '\t'.join(['%s' % len(dict_sVars[sRegion][sClass]) for sClass in list_sOutputKey])
        print('%s\t%s' % (sRegion if sRegion != 'exonic' else 'cds', sOut))
    #loop END: sRegion
#def END: determine_intermediates




def normalize_enrich_byrbp_v1 (sWorkDir, sClipDir, fRecurrency, bPatCnt, bStrictPeaks, bAllRBPs):
    print(green('Determine enrichment by RBP %s' % ('PosxPat' if bPatCnt else 'ByPos'), 1))

    bRecurrency      = 0
    sRBPGroup        = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall        = 'strict' if bStrictPeaks else 'allpeaks'
    sAnalysisTag     =  '%s%s%s%s' % (sSEQDATA, '.recurrent' if bRecurrency else '.nonrecurrent',
                                      '.posxpat' if bPatCnt else '.bypos', '.%s' % sPeakCall)

    dict_sAllVars    = load_all_variant_file (sWorkDir, sRBPGroup, bRecurrency, bPatCnt, sPeakCall)
    dict_nSizeRef    = load_region_sizes     (sWorkDir, sRBPGroup, sPeakCall)

    list_cClipData   = get_clipdata_list(bAllRBPs)
    dict_sRBPID      = {}
    for cClipData in list_cClipData:
        if cClipData.sGeneSym not in dict_sRBPID:
            dict_sRBPID[cClipData.sGeneSym] = []
        dict_sRBPID[cClipData.sGeneSym].append(cClipData.sAccID)
    #loop END: cClipData

    sInDir           = '%s/variant_assigned_byRBP_%s_%s_%s' % (sWorkDir, sPeakCall, sRBPGroup, sSEQDATA)
    sSizeDir         = '%s/%s/%s-RBPsizes'                  % (sClipDir, sPeakCall, sRBPGroup)

    list_sRBPID      = list(dict_sRBPID.keys())
    dict_sOutput     = {}
    list_sRegion     = ['intronic', 'exonic', '5UTR', '3UTR']

    for i, sRBPID in enumerate(list_sRBPID):

        ######################################
        #if sRBPID != 'HNRNPU-human': continue
        ######################################

        sVarFile      = '%s/%s_variant_assigned_byRBP_%s.dat' % (sInDir, sRBPID, sSEQDATA)
        InFile        = open(sVarFile, 'rb')
        dict_sVars    = pickle.load(InFile)
        InFile.close()

        sSizeFile     = '%s/%s_peaksize_combined.dat'         % (sSizeDir, sRBPID)
        InFile        = open(sSizeFile, 'rb')
        dict_nSizeRBP = pickle.load(InFile)
        InFile.close()

        print(sRBPID, i, len(list_sRBPID))
        if sRBPID not in dict_sOutput:
            dict_sOutput[sRBPID] = {}


        dict_sOutput[sRBPID]  = calculate_enrich_by_rbp (list_sRegion, dict_sVars,
                                                         dict_nSizeRBP, dict_nSizeRef, dict_sAllVars,
                                                         bRecurrency, fRecurrency, bPatCnt)
    #loop END: sRBPID
    #output_forPhase2 (sWorkDir, dict_sOutput, bRecurrency, bPatCnt)
    sOutFile     = '%s/%s_%s_fEnrichbyRBP_%s.txt-NumerOnly' % (sWorkDir, sSTUDY, sRBPGroup, sAnalysisTag)
    sHeader      = '\t'.join([sRegion for sRegion in list_sRegion])

    OutFile      = open(sOutFile, 'w')
    OutFile.write('RBPName\t%s\n' % sHeader)
    for sRBP in dict_sOutput:
        sOut = '\t'.join(['%0.2f' % dict_sOutput[sRBP][sRegion] for sRegion in list_sRegion])
        sOut = '%s\t%s\n'         % (sRBP, sOut)
        #print(sOut[:-1])
        OutFile.write(sOut)
    #loop END: sRBP
    OutFile.close()

    sOutFile     = '%s/%s_%s_RBPSizes_%s.txt' % (sWorkDir, sSTUDY, sRBPGroup, sAnalysisTag)
    output_individual_RBP_size(sSizeDir, list_sRBPID, list_sRegion, sOutFile, dict_nSizeRef)

#def END: normalize_enrich_byrbp_old


def qsub_normalize_enrich_byrbp (sWorkDir, sClipDir, fRecurrency, bStrictPeaks, bAllRBPs, sTempScript, sQueue, bTestRun):
    print(green('Determine enrichment by RBP', 1))

    bRecurrency      = 0
    bPatCnt          = 0
    sRBPGroup        = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall        = 'strict' if bStrictPeaks else 'allpeaks'
    sAnalysisTag     =  '%s%s%s%s' % (sSEQDATA, '.recurrent' if bRecurrency else '.nonrecurrent',
                                    '.posxpat' if bPatCnt else '.bypos', '.%s' % sPeakCall)
    sJobName         = 'Jinman.CalcEnrichByRBP.%s.%s.%s'    % (sSTUDY, sPeakCall, sRBPGroup)
    sOutDir          = '%s/variant_assigned_byRBP_%s_%s_%s' % (sWorkDir, sPeakCall, sRBPGroup, sSEQDATA)
    sLogDir          = '%s/log/%s/%s'                       % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sLogDir, exist_ok=True)


    list_cClipData   = get_clipdata_list(bAllRBPs)
    dict_sRBPID      = {}
    for cClipData in list_cClipData:
        if cClipData.sGeneSym not in dict_sRBPID:
            dict_sRBPID[cClipData.sGeneSym] = []
        dict_sRBPID[cClipData.sGeneSym].append(cClipData.sAccID)
    #loop END: cClipData

    sAllVarsFile     = '%s/%s_variant_assigned_%s.data'       % (sWorkDir, sRBPGroup, sAnalysisTag)

    sInDir           = '%s/variant_assigned_byRBP_%s_%s_%s'   % (sWorkDir, sPeakCall, sRBPGroup, sSEQDATA)
    sSizeDir         = '%s/%s/%s-RBPsizes'                    % (sClipDir, sPeakCall, sRBPGroup)
    list_sRBPID      = list(dict_sRBPID.keys())

    for sRBPID in list_sRBPID:
        #InFiles
        sVarFile      = '%s/%s_variant_assigned_byRBP_%s.dat' % (sInDir, sRBPID, sSEQDATA)
        sRBPSizeFile  = '%s/%s_peaksize_combined.dat'         % (sSizeDir, sRBPID)
        sRegionSize   = '%s/%s_RBP_Region_%s_for%s.data'      % (sWorkDir, sRBPGroup, sPeakCall, sSEQDATA)

        #OutFiles
        sOutFile      = '%s/%s_norm_enrich_byRBP_%s.dat'      % (sOutDir, sRBPID, sSEQDATA)
        sSizeOutFile  = '%s/%s_region_size_byRBP_%s.dat'      % (sOutDir, sRBPID, sSEQDATA)

        sScript       = '%s normalize_enrich_byrbp_v2 %s %s %s %s %s %s %s %s %s %s'        \
                         % (sTempScript, sRBPID, sAllVarsFile, sVarFile, sRBPSizeFile, sRegionSize,
                            sOutFile, sSizeOutFile, fRecurrency, bRecurrency, bPatCnt)

        if bTestRun:
            print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s'
                      % (sScript, sLogDir, sQueue, sJobName, sRBPID))
        #if END:
    #loop END: sRBPID
#def END: qsub_normalize_enrich_byrbp


def normalize_enrich_byrbp_v2 (sRBPID, sAllVarFile, sVarFile, sRBPSizeFile, sSizeRegSize, sOutFile, sSizeOutFile, fRecurrency, bRecur, bPatCnt):
    print(green('Determine enrichment by RBP %s' % sRBPID, 1))

    bRecurrency      = int(bRecur)
    fRecurrency      = float(fRecurrency)
    bPatCnt          = int(bPatCnt)

    # Load All VarFile for NoSite
    InFile           = open(sAllVarFile, 'rb')
    dict_sAllVars    = pickle.load(InFile)
    InFile.close()

    # Load RBP VarFile
    InFile           = open(sVarFile, 'rb')
    dict_sVars       = pickle.load(InFile)
    InFile.close()

    # Load RBP SizeFile
    InFile           = open(sRBPSizeFile, 'rb')
    dict_nSizeRBP    = pickle.load(InFile)
    InFile.close()

    # Load Region SizeFile
    InFile           = open(sSizeRegSize, 'rb')
    dict_nSizeRef    = pickle.load(InFile)
    InFile.close()

    list_sRegion     = ['intronic', 'exonic', '5UTR', '3UTR']
    dict_fEnrich     = calculate_enrich_by_rbp (list_sRegion, dict_sVars, dict_nSizeRBP, dict_nSizeRef,
                                                            dict_sAllVars, bRecurrency, fRecurrency, bPatCnt)
    #loop END: sRBPID
    #output_forPhase2 (sWorkDir, dict_sOutput, bRecurrency, bPatCnt)

    sHeader      = '\t'.join([sRegion for sRegion in list_sRegion])
    #OutFile      = open(sOutFile, 'w')
    #OutFile.write('RBPName\t%s\n' % sHeader)
    for sRBP in dict_fEnrich:
        sOut = '\t'.join(['%0.2f' % dict_fEnrich[sRBP][sRegion] for sRegion in list_sRegion])
        sOut = '%s\t%s\n' % (sRBP, sOut)
        print(sOut[:-1])
        #OutFile.write(sOut)
    #loop END: sRBP
    #OutFile.close()

    output_individual_RBP_size (list_sRegion, sSizeOutFile, dict_nSizeRBP, dict_nSizeRef)
#def END: normalize_enrich_byrbp


def calculate_enrich_by_rbp (list_sRegion, dict_sVars, dict_nSizeRBP, dict_nSizeRef, dict_sNoSite, bRecurrency, fRecurrency, bPatCnt):

    list_sOutputKey  = ['RBPOnly', 'MirOnly', 'RBPMir', 'Inter', 'NoSite']

    dict_sLocale     = {sRegion: {sClass:[] for sClass in list_sOutputKey} for sRegion in list_sRegion}

    for sRegion in list_sRegion:
        for sClass in list_sOutputKey:

            dict_sPosKey = dict_sVars[sRegion][sClass]

            if sClass == 'NoSite':
                dict_sPosKey = dict_sNoSite[sRegion][sClass]    

            for nPos in dict_sPosKey:

                nPatCnt = dict_sPosKey[nPos]
                if bRecurrency:
                    if nPatCnt <= int(nSAMPLE_CNT * fRecurrency): continue

                if bPatCnt:
                    for i in range(nPatCnt):
                        dict_sLocale[sRegion][sClass].append(nPos)
                    #loop END: i
                else:
                    dict_sLocale[sRegion][sClass].append(nPos)
                #if END: bPatCnt
            #loop END: nPos
        #loop END: sClass
    #loop END: sRegion

    dict_sOutput = {}
    for sLocale in list_sRegion:

        if sLocale == 'intergenic': continue
        if sLocale == 'other': continue

        if sLocale not in dict_sOutput:
            dict_sOutput[sLocale] = 0.0

        fRefSize             = dict_nSizeRef[sLocale]
        fRBPSize             = dict_nSizeRBP[sLocale]
        dict_sOutput2        = {sClass: 0.0 for sClass in list_sOutputKey}
        for sClass in list_sOutputKey:

            if sClass in ['RBPOnly', 'RBPMir']:
                nDenom = fRBPSize
            else:
                nDenom = fRefSize

            if nDenom == 0:
                print(dict_nSizeRef)
                print(dict_nSizeRBP)
                sys.exit()

            nCategoryCnt = len(dict_sLocale[sLocale][sClass])

            if nCategoryCnt:
                nNormCnt     = (nCategoryCnt / nDenom) * 1000000
                #nNormCnt     = nCategoryCnt
            else: nNormCnt   = 0

            dict_sOutput2[sClass] = nNormCnt
        #loop END: sClass

        try:
            fEnrichRatio = dict_sOutput2['RBPOnly'] / dict_sOutput2['NoSite']
            if fEnrichRatio == 0:
                dict_sOutput[sLocale] = 0
            else:
                #dict_sOutput[sLocale] = math.log(fEnrichRatio,2)
                dict_sOutput[sLocale] = dict_sOutput2['RBPOnly']
        except ZeroDivisionError:
            sys.exit('ZeroDivisionError, %s %s' % (dict_sOutput2['RBPOnly'],dict_sOutput2['NoSite']))
            dict_sOutput[sLocale] = 0
        #try END:
        #sOut     = '\t'.join(['%0.2f' % dict_sOutput2[sClass] for sClass in list_sOutputKey])
        #print('%s\t%s' % (sLocale if sLocale != 'exonic' else 'cds', sOut))
    #loop END:sLocale
    return dict_sOutput
#def END: normalize_enrich_by_rbp


def output_individual_RBP_size (sSizeDir, list_sRBPID, list_sRegion, sOutFile, dict_nSizeRef):

    print(green('Determine Individual RBP Binding Region Sizes', 1))

    dict_fSizes  = {}
    for sRBPID in list_sRBPID:

        if sRBPID not in dict_fSizes:
            dict_fSizes[sRBPID] = {}

        sSizeFile     = '%s/%s_peaksize_combined.dat' % (sSizeDir, sRBPID)
        InFile        = open(sSizeFile, 'rb')
        dict_nSizeRBP = pickle.load(InFile)
        InFile.close()

        for sLocale in list_sRegion:

            if sLocale == 'intergenic': continue
            if sLocale == 'other': continue

            if sLocale not in dict_fSizes[sRBPID]:
                dict_fSizes[sRBPID][sLocale] = ''

            fRefSize = dict_nSizeRef[sLocale]
            fRBPSize = dict_nSizeRBP[sLocale]
            dict_fSizes[sRBPID][sLocale] = '%s-%s-%s' % (fRBPSize, fRefSize-fRBPSize, fRefSize)
    #loop END: sLocale
    #loop END: sRBPID

    sHeader      = '\t'.join([sRegion for sRegion in list_sRegion])
    OutFile      = open(sOutFile, 'w')
    OutFile.write('RBPName\t%s\n' % sHeader)
    for sRBP in dict_fSizes:
        sOut = '\t'.join(['%s' % dict_fSizes[sRBP][sRegion] for sRegion in list_sRegion])
        sOut = '%s\t%s\n' % (sRBP, sOut)
        #print(sOut[:-1])
        OutFile.write(sOut)
    #loop END: sRBP
    OutFile.close()
#def END: output_individual_RBP_size


def load_all_variant_file(sWorkDir, sRBPGroup, bRecurrency, bPatCnt, sPeakCall):
    sInFile = '%s/%s_variant_assigned_%s%s%s%s.data' % (sWorkDir, sRBPGroup, sSEQDATA,
                                                        '.recurrent' if bRecurrency else '.nonrecurrent',
                                                        '.posxpat' if bPatCnt else '.bypos',
                                                        '.%s' % sPeakCall)
    InFile = open(sInFile, 'rb')
    dict_sAllVars = pickle.load(InFile)
    InFile.close()

    return dict_sAllVars
#def END: load_all_variant_file


def load_region_sizes (sWorkDir, sRBPGroup, sPeakCall):
    sRegionSize   = '%s/%s_RBP_Region_%s_for%s.data' % (sWorkDir, sRBPGroup, sPeakCall, sSEQDATA)
    InFile        = open(sRegionSize, 'rb')
    dict_nSizeRef = pickle.load(InFile)
    InFile.close()
    return dict_nSizeRef
#def END: load_region_sizes


def output_enrich_by_rbp (sWorkDir, dict_cVCF, dict_sRBP, dict_sAllVars, dict_nSizeRef, fRecurrency, nBufferSize, nVarWindow):
    print(green('Determine enrichment by RBP', 1))

    bRecurrency  = 0
    bPatCnt      = 0

    sInDir       = '%s/clipdata/byRBPs' % sWorkDir
    list_sRBPs   = list(dict_sRBP.keys())
    list_sNoSite = dict_sAllVars['NoSite']

    dict_sOutput = {}

    for sRBP in list_sRBPs:

        sAssigedFile = '%s/%s_assigned_variants_buffer%s_win%s_top%s.dat' \
                       % (sInDir, sRBP, nBufferSize, nVarWindow, int(fTOP_EXPR * 100))

        sSizeFile    = '%s/%s_rbpsize_buffer%s_win%s_top%s.dat' \
                       % (sInDir, sRBP, nBufferSize, nVarWindow, int(fTOP_EXPR * 100))

        InFile       = open(sAssigedFile, 'rb')
        dict_sVars   = pickle.load(InFile)
        InFile.close()

        InFile        = open(sSizeFile, 'rb')
        dict_nSizeRBP = pickle.load(InFile)
        InFile.close()

        if sRBP not in dict_sOutput:
            dict_sOutput[sRBP] = {}

        dict_sOutput[sRBP] = normalize_enrich_by_rbp(dict_cVCF, dict_sVars, list_sNoSite, dict_nSizeRBP, dict_nSizeRef, bRecurrency, fRecurrency, bPatCnt)
    #loop END: sRBP

    output_forPhase2 (sWorkDir, dict_sOutput, bRecurrency, bPatCnt)


    list_sRegion = ['intronic', 'exonic', '5UTR', '3UTR']

    sOutFile     = '%s/%s_fEnrichbyRBP%s%s.txt' % (sWorkDir, sCELLS,
                                              '.recurrent' if bRecurrency else '.nonrecurrent',
                                              '.patcnt' if bPatCnt else '.single')
    sHeader      = '\t'.join([sRegion for sRegion in list_sRegion])

    OutFile      = open(sOutFile, 'w')
    OutFile.write('RBPName\t%s\n' % sHeader)
    for sRBP in dict_sOutput:
        sOut = '\t'.join(['%0.2f' % dict_sOutput[sRBP][sRegion] for sRegion in list_sRegion])
        sOut = '%s\t%s\n' % (sRBP, sOut)
        #print(sOut[:-1])
        OutFile.write(sOut)
    #loop END: sRBP
    OutFile.close()
#def END: output_enrich_by_rbp


def normalize_enrich (dict_cVCF, dict_sVars, fRecurrency):
    print(green('Determine enrichment', 1))

    bRecurrency      = 0
    bPatCnt          = 0

    list_sOutputKey  = ['RBPOnly', 'MirOnly', 'RBPMir', 'Inter', 'NoSite']
    list_sRegion     = ['intergenic', 'intronic', 'exonic', '5UTR', '3UTR', 'other']

    dict_sLocale     = {sRegion: {sClass:[] for sClass in list_sOutputKey} for sRegion in list_sRegion}

    sHeader          = '\t'.join([sClass for sClass in list_sOutputKey])

    print('region\t%s' % sHeader)
    for sClass in dict_sVars:

        for sPosKey in dict_sVars[sClass]:
            sChrID        = sPosKey.split(',')[0]
            nPos          = int(sPosKey.split(',')[1])
            nPatCnt       = len(dict_cVCF[sChrID][nPos])

            if bRecurrency:
                if nPatCnt <= int(nSAMPLE_CNT * fRecurrency): continue

            list_sLocale  = []
            for cVCF in dict_cVCF[sChrID][nPos]:
                #dict_sVCFInfo = dict(sInfo.split('=') for sInfo in cVCF.sInfo.split(';')[:-1] if len(sInfo.split('=')) == 2)
                #sLocale       = dict_sVCFInfo['Func.refGene'].split('\\x3b')[0]
                sKey          = dict_sGENIC_alt[cVCF.sLocale]
                list_sLocale.append(sKey)
            #loop END:
            sLocale = list(set(list_sLocale))[0]

            #if sLocale == '3UTR' and sClass == 'RBPOnly':
                #print(sClass, sChrID, nPos, nPatCnt)

            if bPatCnt:
                for i in range(len(dict_cVCF[sChrID][nPos])):
                    dict_sLocale[sLocale][sClass].append(nPos)
                #loop END: i

            else:
                dict_sLocale[sLocale][sClass].append(nPos)
            #if END: bPatCnt

        #loop END: sPosKey
    #loop END: sClass

    for sLocale in list_sRegion:
        sOut    = '\t'.join(['%s' % len(dict_sLocale[sLocale][sClass]) for sClass in list_sOutputKey])
        print('%s\t%s' % (sLocale if sLocale != 'exonic' else 'cds', sOut))
    #loop END:sLocale
#def END: normalize_enrich


def normalize_enrich_by_rbp (dict_cVCF, dict_sVars, list_sNoSite, dict_nSizeRBP, dict_nSizeRef, bRecurrency, fRecurrency, bPatCnt):

    list_sOutputKey  = ['RBPOnly', 'MirOnly', 'RBPMir', 'Inter', 'NoSite']
    list_sRegion     = ['intergenic', 'intronic', 'exonic', '5UTR', '3UTR', 'other']

    dict_sLocale     = {sRegion: {sClass:[] for sClass in list_sOutputKey} for sRegion in list_sRegion}

    sHeader          = '\t'.join([sClass for sClass in list_sOutputKey])

    #print('region\t%s' % sHeader)

    dict_sOutput      = {}

    for sClass in dict_sVars:

        list_sPosKey = dict_sVars[sClass]
        if sClass == 'NoSite':
            list_sPosKey = list_sNoSite

        for sPosKey in list_sPosKey:
            sChrID        = sPosKey.split(',')[0]
            nPos          = int(sPosKey.split(',')[1])
            nPatCnt       = len(dict_cVCF[sChrID][nPos])

            if bRecurrency:
                if nPatCnt <= int(nSAMPLE_CNT * fRecurrency): continue

            list_sLocale  = []
            for cVCF in dict_cVCF[sChrID][nPos]:
                sKey          = dict_sGENIC_alt[cVCF.sLocale]
                list_sLocale.append(sKey)
            #loop END:
            sLocale = list(set(list_sLocale))[0]


            if bPatCnt:
                for i in range(len(dict_cVCF[sChrID][nPos])):
                    dict_sLocale[sLocale][sClass].append(nPos)
                #loop END: i

            else:
                dict_sLocale[sLocale][sClass].append(nPos)
            #if END: bPatCnt
        #  #loop END: sPosKey
    #loop END: sClass

    for sLocale in list_sRegion:

        if sLocale == 'intergenic': continue
        if sLocale == 'other': continue

        if sLocale not in dict_sOutput:
            dict_sOutput[sLocale] = 0.0

        sRefSize     = dict_nSizeRef[sLocale]
        sRBPSize     = dict_nSizeRBP[sLocale]

        dict_sOutput2 = {sClass: 0.0 for sClass in list_sOutputKey}

        for sClass in list_sOutputKey:

            if sClass in ['RBPOnly', 'RBPMir']:
                nDenom = sRBPSize
            else:
                nDenom = sRefSize

            nCategoryCnt = len(dict_sLocale[sLocale][sClass])

            if nCategoryCnt:
                nNormCnt     = (nCategoryCnt / nDenom) * 1000000
            else: nNormCnt   = 0

            dict_sOutput2[sClass] = nNormCnt
        #loop END: sClass

        try:
            fEnrichRatio = dict_sOutput2['RBPOnly'] / dict_sOutput2['NoSite']
            if fEnrichRatio == 0:
                dict_sOutput[sLocale] = 0
            else:
                dict_sOutput[sLocale] = math.log(fEnrichRatio,2)
        except ZeroDivisionError:
            dict_sOutput[sLocale] = 0

        sOut     = '\t'.join(['%0.2f' % dict_sOutput2[sClass] for sClass in list_sOutputKey])
        print('%s\t%s' % (sLocale if sLocale != 'exonic' else 'cds', sOut))

    #loop END:sLocale
    return dict_sOutput
#def END: normalize_enrich_by_rbp


def output_forPhase2 (sWorkDir, list_sRegion, dict_fEnrich, bRecurrency, bPatCnt):
    print(green('Output RBP list for Phase 2', 1))

    sOutDir      = '%s/forPhase2' % sWorkDir
    os.makedirs(sOutDir, exist_ok=True)


    for sGenicRegion in list_sRegion:

        sOutFile     = '%s/%s_RBPs_%s%s%s.txt' % (sOutDir, sCELLS, sGenicRegion if sGenicRegion != 'exonic' else 'cds',
                                           '.recurrent' if bRecurrency else '.nonrecurrent',
                                           '.patcnt' if bPatCnt else '.single')

        list_sOutput = [[sRBP, dict_fEnrich[sRBP][sGenicRegion]] for sRBP in dict_fEnrich]
        list_sOutput = sorted(list_sOutput, key=lambda e:e[1], reverse=True)

        OutFile      = open(sOutFile, 'w')
        for sRBP, fEnrich in list_sOutput:
            sOut = '%s\t%s\n' % (sRBP, fEnrich)
            OutFile.write(sOut)
        #loop END: sRBP, fEnrich
        OutFile.close()
    #loop END: sGenicRegion
    sOutFile = '%s/%s_RBP_wEnrich%s%s.dat' % (sWorkDir, sCELLS, '.recurrent' if bRecurrency else '.nonrecurrent',
                                           '.patcnt' if bPatCnt else '.single')
    OutFile  = open(sOutFile, 'wb')
    pickle.dump(dict_fEnrich, OutFile)
    OutFile.close()
#def END: output_forPhase2

## endregion

## endregion

## region RBP-TCGA Phase 2

## region RBPxRBP, Cancer by Cancer analysis

def rbp_by_rbp_analysis (nVarWindow, nClustDist, nBufferSize, fRecurrency, bStrictPeaks, bAllRBPs):

    bPatCnt         = 1
    bRecurrency     = 0
    sRBPGroup       = 'AllRBPs' if bAllRBPs else sCELLS
    sPeakCall       = 'strict' if bStrictPeaks else 'allpeaks'
    sWorkDir        = '%s/03_determine_enrich'       % sBASE_DIR
    sEnrichDir      = '%s/EnrichmentbyRBP'           % sWorkDir
    sOutDir         = '%s/outbyregion_%s_PosxPat_Numer'  % (sEnrichDir, sRBPGroup)
    os.makedirs(sOutDir, exist_ok=True)

    dict_sCancers   = {'TCGA_LIHC': 'Liver',     'TCGA_LAML': 'BoneMarrow',
                       'TCGA_BRCA': 'Breast',    'TCGA_LUAD': 'Lung',
                       'TCGA_UCEC': 'Uterus',    'TCGA_OV'  : 'Ovary',
                       'TCGA_COAD': 'Colorectal','TCGA_CESC': 'Cervix',
                       'TCGA_SARC': 'SoftTissue','TCGA_ESCA': 'Esophagus'}

    list_sRegion    = ['intronic', 'cds', '5UTR', '3UTR']
    dict_sRBPID     = get_clipdata_list_v2(bAllRBPs)

    dict_fEnrich    = load_enrichment_values (sEnrichDir, dict_sRBPID, list_sRegion, bRecurrency,
                                              bPatCnt, sPeakCall, sRBPGroup)

    list_sRBPIDs    = list(dict_sRBPID[sSTUDY])
    list_sTissue    = [dict_sCancers[sCancer] for sCancer in list_sCANCERs]

    for sRegion in list_sRegion:
        dict_sOutput = {sRBPID: {sCancer: dict_fEnrich[sRegion][sCancer][sRBPID] for sCancer in list_sCANCERs} for sRBPID in list_sRBPIDs}
        
        sOutFile = '%s/fEnrichment_%s_%s.txt' % (sOutDir, sRegion, sRBPGroup)
        print('OutFile', sOutFile)

        OutFile  = open(sOutFile, 'w')
        sHeader  = 'CancerTypes\t%s\n' % ('\t'.join(list_sTissue))
        #print(sHeader)
        OutFile.write('%s' % sHeader)
        
        for sRBPID in list_sRBPIDs:
            sOut = '%s\t%s\n' % (sRBPID.replace('-human',''), '\t'.join([dict_sOutput[sRBPID][sCancer] for sCancer in list_sCANCERs]))
            #print(sOut[:-1])
            OutFile.write(sOut)
        #loop END: sRBPID
        OutFile.close()
    #loop END: sRegion
#def END: rbp_by_rbp_analysis


def load_enrichment_values (sEnrichDir, dict_sRBPID, list_sRegion, bRecurrency, bPatCnt, sPeakCall, sRBPGroup):

    dict_sOutput = {sRegion: {sCancer: {sRBPID: '' for sRBPID in  dict_sRBPID[sCancer]}
                              for sCancer in dict_sRBPID} for sRegion in list_sRegion}

    for sCancer in list_sCANCERs:

        if sRBPGroup != 'AllRBPs':
            sCells = 'K562' if sCancer == 'TCGA_LAML' else 'Both'
        else:
            sCells = sRBPGroup

        sInFile      = '%s/%s_%s_fEnrichbyRBP_%s%s%s%s.txt-NumerOnly' % (sEnrichDir, sCancer, sCells, sSEQDATA,
                                                           '.recurrent' if bRecurrency else '.nonrecurrent',
                                                           '.posxpat' if bPatCnt else '.bypos',
                                                           '.%s' % sPeakCall)
        if os.path.isfile(sInFile):
            print('File Found %s' % sInFile)
        else:
            print('****************File Not Found %s' % sInFile)
            sys.exit()
        #if END:

        dict_fEnrich = parse_enrich (sInFile)
        
        for sRBPID in dict_fEnrich:

            for sRegion in dict_fEnrich[sRBPID]:

                fEnrichValue = dict_fEnrich[sRBPID][sRegion]

                dict_sOutput[sRegion][sCancer][sRBPID] = fEnrichValue
            #loop END: sRegion
        #loop END: sRBPID
    #loop END: sCancer
    '''
    for sRegion in dict_sOutput:
        if sRegion != 'cds': continue
        print('sRegion', sRegion)
        for sCancer in dict_sOutput[sRegion]:
            print('sCancer', sCancer)
            if sCancer != 'TCGA_BRCA': continue

            for sRBPID in dict_sOutput[sRegion][sCancer]:

                print(sRBPID, dict_sOutput[sRegion][sCancer][sRBPID])
    '''
    return dict_sOutput
#def END: load_enrichment_values


def parse_enrich (sInFile):

    dict_sOutput = {}
    InFile       = open(sInFile, 'r')

    for sReadLine in InFile:
        '''
        0   sRBPID           RPS11-human                  
        1   fIntron          0.85
        2   fCDS             -0.26
        3   f5UTR            0.41
        4   f3UTR            -0.06
        '''
        if sReadLine.startswith('RBPName'):continue
        list_sColumn  = sReadLine.strip('\n').split('\t')
        sRBPID        = list_sColumn[0]
        fIntron       = list_sColumn[1]
        fCDS          = list_sColumn[2]
        f5UTR         = list_sColumn[3]
        f3UTR         = list_sColumn[4]

        if sRBPID not in dict_sOutput:
            dict_sOutput[sRBPID] = ''
        dict_sOutput[sRBPID] = {'intronic':fIntron, 'cds':fCDS, '5UTR':f5UTR, '3UTR':f3UTR}
    #loop END: sReadLine
    InFile.close()

    #V-S Check
    if not dict_sOutput:
        sys.exit('Invalid Dict : parse_enrich : dict_sOutput size= %d' % len(dict_sOutput))

    return dict_sOutput
#def END: parse_enrich

## endregion

## region Top RBP Examination

def top_RBP_examination (sStudy, nVarWindow, nClustDist, nBufferSize, fRecurrency, bStrictPeaks, bAllRBPs):
    #print(green('Examine Top RBPs', 1))
    bRecurrency      = 0
    bPatCnt          = 1

    sCells           = dict_sCELLS[sStudy]
    sRBPGroup        = 'AllRBPs' if bAllRBPs else sCells
    sPeakCall        = 'strict' if bStrictPeaks else 'allpeaks'

    sWorkDir         = '%s/03_determine_enrich'       % sBASE_DIR
    sStudyDir        = '%s/%s'                        % (sWorkDir, sStudy)
    sClipDir         = '%s/clipdata'                  % sWorkDir
    sPeakDir         = '%s/%s/%s-RBPs_bychr'          % (sClipDir, sPeakCall, sRBPGroup)

    '''
    sInFile          = '%s/%s_variant_assigned_%s%s%s%s.data' % (sWorkDir, sRBPGroup, sSEQDATA,
                                                                '.recurrent' if bRecurrency else '.nonrecurrent',
                                                                '.posxpat' if bPatCnt else '.bypos',
                                                                '.%s' % sPeakCall)
    InFile           = open(sInFile, 'rb')
    dict_sAllVars    = pickle.load(InFile)
    InFile.close()
    '''
    list_cClipData   = get_clipdata_list(bAllRBPs)
    dict_sRBPID      = {}
    for cClipData in list_cClipData:
        if cClipData.sGeneSym not in dict_sRBPID:
            dict_sRBPID[cClipData.sGeneSym] = []
        dict_sRBPID[cClipData.sGeneSym].append(cClipData.sAccID)
    #loop END: cClipData

    sInDir        = '%s/variant_assigned_byRBP_%s_%s_%s' % (sStudyDir, sPeakCall, sRBPGroup, sSEQDATA)
    list_sRBPID   = list(dict_sRBPID.keys())

    ### Manual Inputs ################
    sTargetRBP    = 'EFTUD2-human'
    sTargetRegion = '3UTR'
    #################################

    list_nVarPos          = get_top_candidate_positions (sInDir, sTargetRBP, sTargetRegion)
    sTopWin, dict_sVarPos = determine_RBP_window (sPeakDir, sTargetRBP, list_nVarPos)

    print(sStudy, sTargetRBP, sTargetRegion)
    print('Total VarCnt', len(list_nVarPos))
    print('Highest VarCnt', sTopWin.replace(',','-'), len(dict_sVarPos))
    print('AvgRecurrence', np.mean([dict_sVarPos[nPos] for nPos in dict_sVarPos]))
    for nPos in dict_sVarPos:
        print(nPos, '%s/%s' % (dict_sVarPos[nPos], dict_sSAMPLECNT[sStudy]))
    list_sOverlapPos      = determine_overlapping_RBPs (sPeakDir, sTargetRBP, sTopWin, list_sRBPID, dict_sVarPos)
    visualize_overlaps (sTargetRBP, sTopWin, list_sOverlapPos, dict_sVarPos)

#def END: top_RBP_examination


def get_top_candidate_positions (sInDir, sTargetRBP, sTargetRegion):

    sVarFile      = '%s/%s_variant_assigned_byRBP_%s-RBPOnly.dat' % (sInDir, sTargetRBP, sSEQDATA)
    if not os.path.isfile (sVarFile):
            sys.exit('File Not Found %s' % sVarFile)

    InFile        = open(sVarFile, 'rb')
    dict_sVars    = pickle.load(InFile)
    InFile.close()

    list_nPos = []

    for sChrID in dict_sVars[sTargetRegion]:
        list_nPos += [[sChrID, nPos, dict_sVars[sTargetRegion][sChrID][nPos]]
                      for nPos in dict_sVars[sTargetRegion][sChrID]]
    #loop END: sChrID
    if not list_nPos: sys.exit('Empty List: list_nPos', list_nPos)

    #Top Position by Recurrency (PatCnt)
    return sorted(list_nPos, key=lambda e:(e[2], e[1]), reverse=True)
#def END: get_top_candidate_positions


def determine_RBP_window (sPeakDir, sTargetRBP, list_nVarPos):

    dict_sPosMatch = {}

    for sChrID, nVarPos, nPatCnt in list_nVarPos:

        sPeakFile     = '%s/%s/peaks_%s_%s.dat'    % (sPeakDir, sTargetRBP, sSEQDATA, sChrID)
        if not os.path.isfile (sPeakFile):
            sys.exit('File Not Found %s' % sPeakFile)

        InFile        = open(sPeakFile, 'rb')
        dict_sPeakPos = pickle.load(InFile)
        InFile.close()

        for sPosKey in dict_sPeakPos[sChrID]:
            nRBPStart, nRBPEnd = [int(sPos) for sPos in sPosKey.split(',')]

            sNewPosKey = '%s:%s' % (sChrID, sPosKey)
            if nRBPStart <= nVarPos <= nRBPEnd:
                if sNewPosKey not in dict_sPosMatch:
                    dict_sPosMatch[sNewPosKey] = {}

                if nVarPos not in dict_sPosMatch[sNewPosKey]:
                    dict_sPosMatch[sNewPosKey][nVarPos] = ''
                dict_sPosMatch[sNewPosKey][nVarPos] = nPatCnt
            #if END:
        #loop END: sPosKey
        #print(sPosKey, dict_sPeakPos[sChrID][sPosKey])
    #loop END: sChrID, nVarPos, nPatCnt

    #Find RBP with most variant positions
    list_sWindows   = [[sWindow, len(dict_sPosMatch[sWindow])] for sWindow in dict_sPosMatch]
    list_sTopWindow = sorted(list_sWindows, key=lambda e:e[1], reverse=True)
    sTopRBPPos      = list_sTopWindow[0][0]
    dict_sTopVarPos = dict_sPosMatch[sTopRBPPos]
    return sTopRBPPos, dict_sTopVarPos
#def END: determine_RBP_window


def determine_overlapping_RBPs (sPeakDir, sTargetRBP, sTopWin, list_sRBPID, dict_sVarPos):

    sChrID       = sTopWin.split(':')[0]
    nWinS, nWinE = [int(sPos) for sPos in sTopWin.split(':')[1].split(',')]

    list_sOutput = []
    for sRBPID in list_sRBPID:
        if sRBPID == sTargetRBP: continue

        sPeakFile     = '%s/%s/peaks_%s_%s.dat'    % (sPeakDir, sRBPID, sSEQDATA, sChrID)
        if not os.path.isfile (sPeakFile):
            sys.exit('File Not Found %s' % sPeakFile)

        InFile        = open(sPeakFile, 'rb')
        dict_sPeakPos = pickle.load(InFile)
        InFile.close()

        for sPosKey in dict_sPeakPos[sChrID]:
            nRBPStart, nRBPEnd = [int(sPos) for sPos in sPosKey.split(',')]
            if nRBPStart >= nWinE: continue
            if nRBPEnd   <= nWinS: continue

            list_sOutput.append([sRBPID, nRBPStart, nRBPEnd])
        #loop END: sPosKey
    #loop END: sRBPID

    if not list_sOutput:
        print('No Other Overlapping RBPs')
    else:
        for sRBPID, nRBPStart, nRBPEnd in list_sOutput:
            print('Overlap Found', sRBPID, nRBPStart, nRBPEnd)
    return list_sOutput
#def END: determine_overlapping_RBPs


def visualize_overlaps (sTargetRBP, sTopWin, list_sOverlapPos, dict_sVarPos):

    sChrID         = sTopWin.split(':')[0]
    nWinS, nWinE   = [int(sPos) for sPos in sTopWin.split(':')[1].split(',')]

    list_nRBPStart = [nRBPStart for sRBPID, nRBPStart, nRBPEnd in list_sOverlapPos] + [nWinS]
    list_nRBPEnd   = [nRBPEnd for sRBPID, nRBPStart, nRBPEnd in list_sOverlapPos]   + [nWinE]

    nWinMax        = max(list_nRBPEnd)
    nWinMin        = min(list_nRBPStart)

    nIndexS        = 0
    nIndexE        = nWinMax - nWinMin + 1

    print(nWinMin, nWinMax)
    print(nIndexS, nIndexE)

    list_nVarIndex = [nVarPos - nWinMin for nVarPos in dict_sVarPos]
    list_nRBPIndex = [nRBPPos - nWinMin for nRBPPos in range(nWinS, nWinE+1)]

    list_sRefSeq   = ['-' for i in range(nIndexE)]
    list_sVarSeq   = [' ' for i in range(nIndexE)]

    for nVarIndex in list_nVarIndex:
        list_sVarSeq[nVarIndex] = '*'

    for nRBPIndex in list_nRBPIndex:
        list_sRefSeq[nRBPIndex] = '0'

    print('\'', ''.join(list_sVarSeq))
    print('\'', ''.join(list_sRefSeq), ' %s' % sTargetRBP.replace('-human', ''))

    for sRBPID, nRBPStart, nRBPEnd in list_sOverlapPos:
        list_sRBPSeq    = ['-' for i in range(nIndexE)]
        list_nOverIndex = [nRBPPos - nWinMin for nRBPPos in range(nRBPStart, nRBPEnd)]
        for nOverIndex in list_nOverIndex:
            list_sRBPSeq[nOverIndex] = '0'
        print('\'', ''.join(list_sRBPSeq), ' %s' % sRBPID.replace('-human', ''))


#def END: visualize_overlaps
## endregion




## region Determine Target Genes
def determine_target_genes (nClustDist, nBufferSize, fRecurrency):
    print(blue('Determine target genes', 1))

    nVarWindow   = 0
    bRecurrency  = 0
    bPatCnt      = 0
    nTopRBPs     = 20  #Top 10 or 20 RBPs by enrichment in region of interest (CDS, 3UTR, or Introns)

    sPrev2Dir    = '%s/02_capture_filter/%s'   % (sBASE_DIR, sSTUDY)
    sPrevDir     = '%s/03_determine_enrich/%s' % (sBASE_DIR, sSTUDY)
    sWorkDir     = '%s/04_targetprofile/%s'    % (sBASE_DIR, sSTUDY)
    os.makedirs(sWorkDir, exist_ok=True)

    #Load RefSeq data
    list_cRefSeq  = load_refseq_data (sPrev2Dir)
    print('list_cRefSeq', len(list_cRefSeq))

    #Acquire top expressed genes
    dict_sGeneSym = acquire_top_expressed_genes (list_cRefSeq)
    print('dict_sGeneSym', len(dict_sGeneSym))

    #Load RBP lists
    dict_sRBP     = load_rbp_lists (sPrevDir, bRecurrency, bPatCnt, nTopRBPs)
    print('dict_sRBP', len(dict_sRBP))

    #Load VCF data
    dict_cVCF     = load_vcf_data  (sPrevDir)
    print('dict_cVCF', len(dict_cVCF))

    #Load Peak Data and determine targets
    dict_sRegion  = load_peak_find_targets_data (sWorkDir, sPrevDir, dict_sRBP, dict_cVCF, nTopRBPs, nBufferSize, nVarWindow)
    print('dict_sRegion', len(dict_sRegion))

    #Profile targets
    profile_rbp_targets (sWorkDir, sPrevDir, dict_sRegion, nTopRBPs, nBufferSize, nVarWindow)
#def END: determine_target_genes


def load_rbp_lists (sPrevDir, bRecurrency, bPatCnt, nTopRBPs):
    print(green('Loading RBP data with enrichment values', 1))

    list_sRegion = ['intronic', 'exonic', '5UTR', '3UTR']
    sInFile      = '%s/%s_RBP_wEnrich%s%s.dat' % (sPrevDir, sCELLS, '.recurrent' if bRecurrency else '.nonrecurrent',
                                           '.patcnt' if bPatCnt else '.single')
    InFile       = open(sInFile, 'rb')
    dict_fEnrich = pickle.load(InFile)
    InFile.close()

    dict_sOutput = {}
    for sGenicRegion in list_sRegion:
        if sGenicRegion not in dict_sOutput:
            dict_sOutput[sGenicRegion] = []

        list_sOutput = [[sRBP, dict_fEnrich[sRBP][sGenicRegion]] for sRBP in dict_fEnrich]
        list_sOutput = sorted(list_sOutput, key=lambda e:e[1], reverse=True)[:nTopRBPs]

        dict_sOutput[sGenicRegion] = list_sOutput
    #loop END: sGenicRegion

    #VS-Check
    if not dict_sOutput:
        sys.exit('Invalid Dict : load_rbp_lists : dict_sOutput size= %d' % len(dict_sOutput))

    return dict_sOutput
#def END: load_rbp_lists


def load_vcf_data (sPrevDir):
    print(green('Loading VCF data', 1))

    sInFile = '%s/VCFs-%s/vcfs_capture.dat'  % (sPrevDir, sSEQDATA)
    InFile  = open(sInFile, 'rb')
    dict_sOutput = pickle.load(InFile)
    InFile.close()

    #VS-Check
    if not dict_sOutput:
        sys.exit('Invalid Dict : load_vcf_data : dict_sOutput size= %d' % len(dict_sOutput))

    return dict_sOutput
#def END: load_vcf_data


def load_peak_find_targets_data (sWorkDir, sPrevDir, dict_sRBP, dict_cVCF, nTopRBPs, nBufferSize, nVarWindow):
    print(green('Load peaks and determine rbp targets', 1))

    sInPeakDir       = '%s/clipdata/byRBPs' % sPrevDir
    list_sRegion     = ['exonic', 'intronic', '5UTR', '3UTR']
    sRegionFile      = '%s/targetgenes_rbps%s_buffer%s_win%s_top%s.dat' % \
                       (sWorkDir, nTopRBPs, nBufferSize, nVarWindow, fTOP_EXPR)

    if not os.path.isfile(sRegionFile):

        dict_sRegion     = {}
        for sRegion in list_sRegion:

            if sRegion not in dict_sRegion:
                dict_sRegion[sRegion] = ''

            dict_sRBPTargets = {}
            dict_sRBPTargets2 = {}

            for sRBP, fEnrich in dict_sRBP[sRegion]:

                if sRBP not in dict_sRBPTargets:
                    dict_sRBPTargets[sRBP] = ''

                sPeakFile    = '%s/%s_assigned_variants_buffer%s_win%s_top%s.dat' \
                                  % (sInPeakDir, sRBP, nBufferSize, nVarWindow, int(fTOP_EXPR*100))
                #V-S Check
                if not os.path.isfile(sPeakFile):
                    sys.exit('File Not Found', sPeakFile)

                InFile       = open(sPeakFile, 'rb')
                dict_sVars   = pickle.load(InFile)
                InFile.close()

                dict_sTargetGenes = {}
                for sPosKey in dict_sVars['RBPOnly']:
                    sChrID        = sPosKey.split(',')[0]
                    nPos          = int(sPosKey.split(',')[1])

                    list_sLocale  = []
                    list_sGeneSym = []
                    for cVCF in dict_cVCF[sChrID][nPos]:
                        list_sLocale.append(dict_sGENIC_alt[cVCF.sLocale])
                        list_sGeneSym.append(cVCF.sGeneSym)
                    #loop END:

                    sLocale     = list(set(list_sLocale))[0]
                    sGeneSym    = list(set(list_sGeneSym))[0].upper()

                    if sLocale != sRegion: continue

                    if sGeneSym not in dict_sRBPTargets2:
                        dict_sRBPTargets2[sGeneSym] = []
                    dict_sRBPTargets2[sGeneSym] += dict_cVCF[sChrID][nPos]

                    if sGeneSym not in dict_sTargetGenes:
                        dict_sTargetGenes[sGeneSym] = []
                    dict_sTargetGenes[sGeneSym] += dict_cVCF[sChrID][nPos]
                #loop END: sPosKey

                dict_sRBPTargets[sRBP] = dict_sTargetGenes

                #print('%s Target Genes %s' % (sRBP, len(dict_sRBPTargets[sRBP])))
            #loop END: sRBP, fEnrich

            print('%s Target Genes %s' % (sRegion, len(dict_sRBPTargets2)))

            dict_sRegion[sRegion] = dict_sRBPTargets
        #loop END: sRegion

        OutFile = open(sRegionFile, 'wb')
        pickle.dump(dict_sRegion, OutFile)
        OutFile.close()
        return dict_sRegion
    else:
        InFile = open(sRegionFile, 'rb')
        dict_sRegion = pickle.load(InFile)
        InFile.close()

        return dict_sRegion
    #if END:
#def END: load_peak_find_targets_data


def profile_rbp_targets (sWorkDir, sPrevDir, dict_sRegion, nTopRBPs, nBufferSize, nVarWindow):
    print(green('Profile RBP targets', 1))

    nVarLimit    = 1

    sOutDir      = '%s/forGOAnalysis'  % sWorkDir
    os.makedirs(sOutDir, exist_ok=True)

    list_sRegion = ['exonic', 'intronic', '5UTR', '3UTR']
    list_sRegion = ['3UTR']
    for sRegion in list_sRegion:

        sOutFile = '%s/%s_rbps%s_buffer%s_win%s_top%s_var%s.txt' % (sOutDir, sRegion, nTopRBPs, nBufferSize, nVarWindow, fTOP_EXPR, nVarLimit)

        dict_sTargetGenes = {}
        dict_sTopTarget   = {}

        for sRBP in dict_sRegion[sRegion]:

            for sTargetGene in dict_sRegion[sRegion][sRBP]:

                if len(dict_sRegion[sRegion][sRBP][sTargetGene]) <= nVarLimit: continue

                #print(sRBP, sTargetGene, len(dict_sRegion[sRegion][sRBP][sTargetGene]))

                list_cVCF   = dict_sRegion[sRegion][sRBP][sTargetGene]

                if sTargetGene not in dict_sTargetGenes:
                    dict_sTargetGenes[sTargetGene] = []
                dict_sTargetGenes[sTargetGene] += list_cVCF

                if sTargetGene not in dict_sTopTarget:
                    dict_sTopTarget[sTargetGene] = [[], []]
                dict_sTopTarget[sTargetGene][0] += list_cVCF
                dict_sTopTarget[sTargetGene][1].append(sRBP)

            #loop END: sTargetGene
        #loop END: sRBP

        list_nVarCnt   = [len(dict_sTargetGenes[sTargetGene]) for sTargetGene in dict_sTargetGenes]
        fMin           = min(list_nVarCnt)
        fMax           = max(list_nVarCnt)
        fMean          = np.mean(list_nVarCnt)
        fMedian        = np.median(list_nVarCnt)


        print('Target Genes in %s by top %s RBPs with >%s:%s,min=%d,max=%d,mean=%d,median=%d'
              % (sRegion, nTopRBPs, nVarLimit, len(dict_sTargetGenes), fMin, fMax, fMean, fMedian))

        list_sTopTarget = [[sTargetGene, len(dict_sTopTarget[sTargetGene][0]), dict_sTopTarget[sTargetGene][1]] for sTargetGene in dict_sTopTarget]
        list_sTopTarget = sorted(list_sTopTarget, key=lambda  e:e[1], reverse=True)

        for sTargetGene, nVarCnt, list_sRBPs in list_sTopTarget:
            print(sTargetGene, nVarCnt, len(dict_sTargetGenes[sTargetGene]))

        sys.exit()

        OutFile = open(sOutFile, 'w')
        for sTargetGene in dict_sTargetGenes:
            sOut = '%s\n' % sTargetGene
            OutFile.write(sOut)
        OutFile.close()
    #loop END: sRegion
#def END: profile_rbp_targets

## endregion


## endregion


def pathway_survey (list_cClipData):

    dict_cClipData = {}
    for cClipData in list_cClipData:

        sKey = cClipData.sGeneSym.split('-')[0]
        if sKey not in dict_cClipData:
            dict_cClipData[sKey] = []
        dict_cClipData[sKey].append(cClipData.sAccID)
    #loop END: cClipData
    print(len(dict_cClipData))

    list_sRBPIDs   = []
    list_sOutput   = []
    dict_sCPDBGene = load_CPDB()
    for sGeneSym in dict_cClipData:

        try: list_sPathways = (dict_sCPDBGene[sGeneSym])
        except KeyError: continue
        # Translation related
        list_sPathways = [sPathway for sPathway in list_sPathways if 'translation' in sPathway or
                          'ribosome' in sPathway]

        if list_sPathways:
            list_sOutput += dict_cClipData[sGeneSym]
            list_sRBPIDs.append(sGeneSym)
            print('%s:%s' % (sGeneSym, list_sPathways[0]))
        #if END:
    #loop END: sGeneSym
    print('Translation Related RBPS: ', len(list_sRBPIDs))

    return list_sOutput
#def END: pathway_survey

def basic_stats_rbp (list_cClipData):
    print(green('PRE-STEP 1: Basic Statistics of RBP Data (eCLIP and iCLIP)', 1))
    bCapture       = 0
    dict_sCapture  = parse_capture_bedfile(sCAPTUREFILE)
    dict_cRefSeq   = cRef_parse_refflat_line(sREFSEQFILE, dict_sCapture, bCapture)

    dict_cRef      = {}
    for sNMID in dict_cRefSeq:

        sKey = '%s,%s' % (dict_cRefSeq[sNMID].sChrID, dict_cRefSeq[sNMID].sStrand)
        if sKey not in dict_cRef:
            dict_cRef[sKey] = []
        dict_cRef[sKey].append(dict_cRefSeq[sNMID])
    #loop END: sNMID

    sClipBedDir    = '%s/sorted_bed'             % sECLIP_DIR
    sOutDir        = '%s/01_basicstats/ClipData' % sBASE_DIR
    os.makedirs(sOutDir, exist_ok=True)

    for cClipData in list_cClipData:

        #if os.path.isfile('%s/%s.peakdata' % (sOutDir, cClipData.sAccID)): continue
        sBedFile        = '%s/%s.bed' % (sClipBedDir, cClipData.sAccID)
        list_cPeakData  = cPeakData_parse_bedfile(sBedFile, dict_cRef)

        print(sBedFile, len(list_cPeakData))

        OutFile         = open('%s/%s.peakdata' % (sOutDir, cClipData.sAccID), 'wb')
        pickle.dump(list_cPeakData, OutFile)
        OutFile.close()
    # loop END: cClipData

    get_peak_distribution (sOutDir, list_cClipData)

    print_done_msg('Basic Stats for RBP')
#def END: basic_stats_rbps


def buffer_rbps (list_cClipData, nBufferSize):
    print(green('PRE-STEP 1-2: Buffering RBPs by BufferSize', 1))
    bCapture       = 0
    dict_sCapture  = parse_capture_bedfile(sCAPTUREFILE)
    dict_cRefSeq   = cRef_parse_refflat_line(sREFSEQFILE, dict_sCapture, bCapture)

    dict_cRef      = {}
    for sNMID in dict_cRefSeq:

        sKey = '%s,%s' % (dict_cRefSeq[sNMID].sChrID, dict_cRefSeq[sNMID].sStrand)
        if sKey not in dict_cRef:
            dict_cRef[sKey] = []
        dict_cRef[sKey].append(dict_cRefSeq[sNMID])
    #loop END: sNMID

    sClipBedDir    = '%s/sorted_bed'                            % sECLIP_DIR
    sOutDir        = '%s/01_basicstats/ClipData_buffered/%sbps' % (sBASE_DIR, nBufferSize)
    os.makedirs(sOutDir, exist_ok=True)

    for cClipData in list_cClipData:

        #if os.path.isfile('%s/%s.peakdata' % (sOutDir, cClipData.sAccID)): continue
        sBedFile        = '%s/%s.bed' % (sClipBedDir, cClipData.sAccID)
        list_cPeakData  = cPeakData_parse_bedfile(sBedFile, dict_cRef, 1)

        print(sBedFile, len(list_cPeakData))
        OutFile         = open('%s/%s.bed' % (sOutDir, cClipData.sAccID), 'w')
        for cPeak in list_cPeakData:
            nNewStart = cPeak.nStartPos - nBufferSize
            nNewEnd   = cPeak.nEndPos   + nBufferSize
            sNewName  = 'Buff-%s'   % cPeak.sSample
            sOut      = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                        (cPeak.sChrID, nNewStart, nNewEnd, sNewName, cPeak.nScore,
                         cPeak.sStrand, cPeak.fSignal, cPeak.fPvalue, '-1', '-1')
            OutFile.write(sOut)
        #loop END: cPeak
        OutFile.close()
    # loop END: cClipData
#def END: buffer_rbps


def basic_stats_vcf_old ():
    print(green('PRE-STEP 2-1: Bastic Stats of VCF data with ANNOVAR annotation information', 1))
    '''Load RefSeq n=18625'''
    bCapture       = 0
    dict_sCapture  = parse_capture_bedfile(sCAPTUREFILE)
    dict_cRefSeq   = cRef_parse_refflat_line(sREFSEQFILE, dict_sCapture, bCapture)

    '''Deteremine regional sizes for normalization '''
    dict_sGenicSize = calculate_genic_size (dict_cRefSeq, dict_sCapture, bCapture)
    sVCFDir         = '%s/01_basicstats/VCFs'      % sBASE_DIR
    sOutDir         = '%s/01_basicstats'           % sBASE_DIR

    # Load Patient IDs
    list_sPatIDs  = get_patientIDs(sSTUDY, sSEQDATA)
    dict_sLocale  = {}
    nTotalVar     = 0

    OutFile       = open('%s/VariantsPerPatients.txt' % (sOutDir), 'w')
    dict_sVarDist = {}
    list_tLOD     = []
    for sPatID in list_sPatIDs:
        sVCFFile  = '%s/%s.hg19_multianno.vcf'  % (sVCFDir, sPatID)
        print(sVCFFile)
        list_cVCF = cVCF_parse_vcf_files (sVCFFile)
        nTotalVar += len(list_cVCF)

        OutFile.write('%s\t%s\n' % (sPatID, len(list_cVCF)))
        if sPatID not in dict_sVarDist:
            dict_sVarDist[sPatID] = {}

        for cVCF in list_cVCF:
            #Split the Info column by ';' then key:value by '='
            dict_sVCFInfo = dict(sInfo.split('=') for sInfo in cVCF.sInfo.split(';')[:-1] if len(sInfo.split('=')) == 2)

            sLocale       = dict_sVCFInfo['Func.refGene'].split('\\x3b')[0]
            sKey          = dict_sGENIC_alt[sLocale]

            if sKey not in dict_sVarDist[sPatID]:
                dict_sVarDist[sPatID][sKey] = 0
            dict_sVarDist[sPatID][sKey] += 1

            if sKey not in dict_sLocale:
                dict_sLocale[sKey] = 0
            dict_sLocale[sKey] += 1

            list_tLOD.append(float(dict_sVCFInfo['TLOD']))
        #loop END: cVCF

    #loop END: sPatID
    OutFile.close()

    nTotalVar2 = 0
    for sLocale in dict_sLocale:
        nTotalVar2 += dict_sLocale[sLocale]
        print(sLocale, (dict_sLocale[sLocale]))
        #print(sLocale, (dict_sLocale[sLocale] / dict_sGenicSize[sLocale]) * 1000000)
    sys.exit()

    print(nTotalVar)
    print(nTotalVar2)
    print('Size', len(list_tLOD))
    print('Min TLOD', min(list_tLOD))
    print('MAX TLOD', max(list_tLOD))

    # VarDist by Patient
    list_sLocale = ['exonic', 'UTRs', 'intergenic', 'intronic']
    OutFile2     = open('%s/VariantsDistPerPatients.txt' % (sOutDir), 'w')
    for sPatID in dict_sVarDist:
        list_nVarCnts = [dict_sVarDist[sPatID][sLocale] for sLocale in list_sLocale]
        sOut          = '%s\t%s\n' % (sPatID, '\t'.join([str(nCnt) for nCnt in list_nVarCnts]))
        OutFile2.write(sOut)
    #loop END: sPatID
    OutFile2.close()

    # tLOD statistics
    print('Size', len(list_tLOD))
    print('Min TLOD', min(list_tLOD))
    print('MAX TLOD', max(list_tLOD))
#def END: basic_stats_vcf


def basic_stats_vcf_v2 ():
    print(green('PRE-STEP 2-1: Bastic Stats of VCF data with ANNOVAR annotation information', 1))

    sVCFDir     = '%s/vardata/%s/VCFs_hg19'            % (sBASE_DIR, sSTUDY)
    # Load Patient IDs
    list_sPatIDs  = get_patientIDs(sSTUDY, sSEQDATA)
    list_tLOD     = []
    list_nLOD     = []

    for sPatID in list_sPatIDs:
        sVCFFile  = '%s/%s.vcf'        % (sVCFDir, sPatID)
        print(sVCFFile)
        list_cVCF = cVCF_parse_vcf_files (sVCFFile)

        for cVCF in list_cVCF:
            #Split the Info column by ';' then key:value by '='
            dict_sVCFInfo = dict(sInfo.split('=') for sInfo in cVCF.sInfo.split(';') if len(sInfo.split('=')) == 2)
            list_tLOD.append(float(dict_sVCFInfo['TLOD']))
            list_nLOD.append(float(dict_sVCFInfo['NLOD']))
        #loop END: vVCF
    #loop END: sPatID

    # tLOD statistics
    print('Size TLOD', len(list_tLOD))
    print('Min TLOD', min(list_tLOD))
    print('MAX TLOD', max(list_tLOD))
    print('Size NLOD', len(list_nLOD))
    print('Min NLOD', min(list_nLOD))
    print('MAX NLOD', max(list_nLOD))
#def END: basic_stats_vcf


def basic_stats_clusters (nClustDist):
    print(green('STEP 2-1: Basic statistics of clusters', 1))

    print_start_msg('Clusters Basic Stats')
    sClusterDir     = '%s/02_rbp_vardata'                % sBASE_DIR
    sClusterFile    = '%s/FullVCFs.clustered.%sbps.txt'  % (sClusterDir, nClustDist)
    sOutDir         = '%s/03_hotspots/%sbps'             % (sBASE_DIR, nClustDist)
    os.makedirs(sOutDir, exist_ok=True)

    list_cClusters  = cClusterData_load_clustered_VCFs (sClusterFile)
    nTotalVarCnt    = len(list_cClusters)
    print('Total Variants', nTotalVarCnt)

    ## Group by cluster IDs
    dict_nClusterID = {}
    for cClust in list_cClusters:

        if cClust.nClusterID not in dict_nClusterID:
            dict_nClusterID[cClust.nClusterID] = []

        dict_nClusterID[cClust.nClusterID].append(cClust)
    #loop END: cVCF
    print('All Clusters', len(dict_nClusterID))

    # Filter out clusters with <2 mutations
    dict_nClusterID = {nID: dict_nClusterID[nID] for nID in dict_nClusterID if len(dict_nClusterID[nID]) > nMIN_VARIANTS}

    print('Filtered Clusters', len(dict_nClusterID))

    list_nSize      = []
    list_nVarCnt    = []
    dict_sTest      = {}
    for nClusterID in dict_nClusterID:
        nStartPos   = min([cClust.nPos for cClust in dict_nClusterID[nClusterID]])
        nEndPos     = max([cClust.nPos for cClust in dict_nClusterID[nClusterID]])
        sChrID      = list(set(([cClust.sChrID for cClust in dict_nClusterID[nClusterID]])))[0]
        nSize       = nEndPos - nStartPos + 1

        #if nSize == 1536:                            print('Largest', len(dict_nClusterID[nClusterID]))
        #if len(dict_nClusterID[nClusterID]) == 8152: print('MostVar', nSize)

        if nClusterID == 1113596:

            print(sChrID, nStartPos, nEndPos, nSize, len(dict_nClusterID[nClusterID]))
            list_sPatID  = list(set(([cClust.sPatID for cClust in dict_nClusterID[nClusterID]])))
            print('Patients', len(list_sPatID))


        if len(dict_nClusterID[nClusterID]) not in dict_sTest:
            dict_sTest[len(dict_nClusterID[nClusterID])] = ''
        dict_sTest[len(dict_nClusterID[nClusterID])] = nClusterID

        list_nSize.append(nSize)
        list_nVarCnt.append(len(dict_nClusterID[nClusterID]))
    #loop END: nClusterID
    list_nVarCnt = sorted(list_nVarCnt, reverse=True)
    print('BasicStats:')
    print('MinSize', min(list_nSize))
    print('MaxSize', max(list_nSize))
    print('AverageSize', np.mean(list_nSize))
    print('MedianSize', np.median(list_nSize))
    print('MinVarCnt', min(list_nVarCnt))
    print('MaxVarCnt', max(list_nVarCnt))
    print('AverageVarCnt', np.mean(list_nVarCnt))
    print('MedianVarCnt', np.median(list_nVarCnt))

    print_done_msg('Clusters Basic Stats')
#def END: basic_stats_clusters


def qsub_evaluate_clusters (nClustDist, nBackWinSize, sTempScript, sQueue, bTestRun, nJobs):
    print(green('STEP 3: Statistically evaluate clusters for hotspots', 1))
    sClusterDir     = '%s/02_rbp_vardata/'                      % sBASE_DIR
    sClusterFile    = '%s/FullVCFs.clustered.%sbps.txt'         % (sClusterDir, nClustDist)
    sOutDir         = '%s/03_hotspots/%sbps/%s'                 % (sBASE_DIR, nClustDist, nBackWinSize)
    sJobName        = 'Jinman.Evaluate.Cluster.%sbps.%s'        % (nClustDist, nBackWinSize)
    sLogDir         = '%s/log/%s/%s'                            % (sBASE_DIR, sJobName, sTIME_STAMP)

    os.makedirs(sLogDir, exist_ok=True)
    os.makedirs(sOutDir, exist_ok=True)

    list_cClusters  = cClusterData_load_clustered_VCFs (sClusterFile)

    ## Group by cluster and patient IDs
    dict_nClusterID = {}
    for cClust in list_cClusters:
        if cClust.nClusterID not in dict_nClusterID:
            dict_nClusterID[cClust.nClusterID] = []
        dict_nClusterID[cClust.nClusterID].append(cClust)
    #loop END: cClust

    # Filter out clusters with <2 mutations
    dict_nClusterID = {nID: dict_nClusterID[nID] for nID in dict_nClusterID if len(dict_nClusterID[nID]) > nMIN_VARIANTS}
    #######################################

    OutFile         = open('%s.%s.dat' % (sClusterFile[:-4], nBackWinSize), 'wb')
    pickle.dump(dict_nClusterID, OutFile)
    OutFile.close()

    nTotalClusters  = len(dict_nClusterID)

    for nJobIndex in range(nJobs):
        ## MD
        nStartIndex      = int(nTotalClusters*(nJobIndex+0)/nJobs)
        nEndIndex        = int(nTotalClusters*(nJobIndex+1)/nJobs)  - 1
        sOutFile         = '%s/%04d.dat'                  % (sOutDir, nJobIndex+1)
        sScript          = 'date; %s evaluate_clusters'   % sTempScript
        sScript         += ' %s %s %s %s %s; date'     % (nClustDist, nStartIndex, nEndIndex,
                                                             nBackWinSize, sOutFile)
        sScript          = 'echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%04d'\
                            % (sScript, sLogDir, sQueue, sJobName, nJobIndex+1)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    #loop END: nJobIndex
#def END: qsub_evaluate_clusters


def evaluate_clusters (nClustDist, nStart, nEnd, nBackWinSize, sOutFile):
    sVCFDir         = '%s/vardata/%s/VCFs_hg19'                     % (sBASE_DIR, sSTUDY)
    sClusterDir     = '%s/02_rbp_vardata'                   % sBASE_DIR
    sClusterFile    = '%s/FullVCFs.clustered.%sbps.%s.dat'  % (sClusterDir, nClustDist, nBackWinSize)
    sKeyFile        = '%s/%s/%s_rnaseq.txt'                 % (sBASE_DIR, sSTUDY, sSTUDY)

    InFile          = open(sClusterFile, 'rb')
    dict_nClusterID = pickle.load(InFile)
    InFile.close()
    '''
    list_nVarCnt      = []
    for nClusterID in dict_nClusterID:

        cHotspot                = cHotspotData()
        cHotspot.nClusterID     = nClusterID
        cHotspot.list_cClusters = dict_nClusterID[nClusterID]

        #All variants in the hotspot, same chromosome id
        #print(list(set([cClust.sChrID for cClust in cHotspot.list_cClusters])))
        assert len(list(set([cClust.sChrID for cClust in cHotspot.list_cClusters]))) == 1

        cHotspot.nStartPos     = min([cClust.nPos for cClust in cHotspot.list_cClusters])
        cHotspot.nEndPos       = max([cClust.nPos for cClust in cHotspot.list_cClusters])
        cHotspot.nSize         = cHotspot.nEndPos - cHotspot.nStartPos + 1
        cHotspot.nVarCnt = len(cHotspot.list_cClusters)
        list_nVarCnt.append(cHotspot.nVarCnt / cHotspot.nSize)

    print('FullCnt', len(list_nVarCnt))
    print('Mean', np.mean(list_nVarCnt))
    print('Max', max(list_nVarCnt))
    print('Min', min(list_nVarCnt))
    sys.exit()
    '''
    fAverageVarCnt  = 5.34429032507

    list_nIDs       = sorted(list(dict_nClusterID.keys()))[int(nStart):int(nEnd)]
    list_cHotspots  = []

    for nClusterID in list_nIDs:

        cHotspot                = cHotspotData()
        cHotspot.nClusterID     = nClusterID
        cHotspot.list_cClusters = dict_nClusterID[nClusterID]

        #All variants in the hotspot, same chromosome id
        #print(list(set([cClust.sChrID for cClust in cHotspot.list_cClusters])))
        assert len(list(set([cClust.sChrID for cClust in cHotspot.list_cClusters]))) == 1
        cHotspot.sChrID        = list(set([cClust.sChrID   for cClust in cHotspot.list_cClusters]))[0]
        cHotspot.sLocale       = list(set([cClust.sLocale  for cClust in cHotspot.list_cClusters]))[0]
        cHotspot.sGeneSym      = list(set([cClust.sGeneSym for cClust in cHotspot.list_cClusters]))[0]
        cHotspot.sGeneID       = list(set([cClust.sGeneID for cClust in cHotspot.list_cClusters]))[0]
        cHotspot.list_sPatIDs  = list(set([cClust.sPatID for cClust in cHotspot.list_cClusters]))
        cHotspot.list_nPos     = list(set([cClust.nPos for cClust in cHotspot.list_cClusters]))
        cHotspot.list_sCovered = check_rnaseq_coverage(sKeyFile, cHotspot.list_sPatIDs)
        cHotspot.nStartPos     = min([cClust.nPos for cClust in cHotspot.list_cClusters])
        cHotspot.nEndPos       = max([cClust.nPos for cClust in cHotspot.list_cClusters])
        cHotspot.nSize         = cHotspot.nEndPos - cHotspot.nStartPos + 1
        cHotspot.nPosCnt       = len(cHotspot.list_nPos)
        cHotspot.nVarCnt       = len(cHotspot.list_cClusters)
        cHotspot.nTotVarCnt    = int(np.mean([cClust.nVarCnt for cClust in cHotspot.list_cClusters]))

        if nBackWinSize != 'genome':
            cHotspot.fBackProb = cHotspot_determine_back_probability(sVCFDir, cHotspot.list_cClusters, cHotspot.sChrID,
                                                                    int(nBackWinSize), cHotspot.nStartPos, cHotspot.nEndPos)
        else:
            cHotspot.fBackProb = np.mean([cClust.nVarCnt / nGENOME_SIZE for cClust in cHotspot.list_cClusters])

        cHotspot.fMeanCnt     = np.mean([cClust.nVarCnt / cHotspot.nSize for cClust in cHotspot.list_cClusters])
        cHotspot.fPvalue      = stats.binom_test(int(cHotspot.nVarCnt / cHotspot.nSize), cHotspot.nTotVarCnt, cHotspot.fBackProb)
        cHotspot.fPvalue2     = stats.distributions.poisson.pmf(int(cHotspot.nVarCnt / cHotspot.nSize), cHotspot.fBackProb)

        list_cHotspots.append(cHotspot)
        '''
        print('*********************************************\t%s\t%s\t%s\t%s' % (nClusterID,
                                                                            'REJECT' if cHotspot.fPvalue > 0.05 else '',
                                                                            'REJECT' if cHotspot.fPvalue2 > 0.05 else '',
                                                                            'REJECT' if cHotspot.fPvalue3 > 0.05 else ''))
        print('cHotspot.fPvalue', cHotspot.fPvalue)
        print('cHotspot.fPvalue2', cHotspot.fPvalue2)
        print('cHotspot.fPvalue3', cHotspot.fPvalue3)
        print('Stats', cHotspot.nSize, cHotspot.nVarCnt, cHotspot.nTotVarCnt, cHotspot.fBackProb)
        '''
    #loop END: nClusterID

    list_cHotspots = sorted(list_cHotspots, key=lambda c: c.fPvalue)

    #V-S Check
    if not list_cHotspots:
        sys.exit('Invalid list : list_cHotspots size= %d' % len(list_cHotspots))

    # OutFile for cluster count distribution per sample
    OutFile  = open(sOutFile, 'wb')
    pickle.dump(list_cHotspots, OutFile)
    OutFile.close()
#def END: evalute_clusters


def check_rnaseq_coverage (sInFile, list_sPatIDs):
    InFile       = open(sInFile, 'r')
    list_sOutput = [sReadLine.strip('\n').split('\t')[3] for sReadLine in InFile]
    InFile.close()
    return list(set(list_sPatIDs) & set(list_sOutput))
#def END: check_rnaseq_coverage


def qsub_multipletest (nClustDist, nBackWinSize, sTempScript, sQueue, bTestRun, nJobs):
    print(green('STEP 4: Perform multiple test correction (FDR)', 1))
    sInDir          = '%s/03_hotspots/%sbps/%s'            % (sBASE_DIR, nClustDist, nBackWinSize)
    sOutDir         = '%s/04_mtcorrected/%sbps/%s'         % (sBASE_DIR, nClustDist, nBackWinSize)
    sJobName        = 'Jinman.MT.Correction.%sbps.%s'      % (nClustDist, nBackWinSize)
    sLogDir         = '%s/log/%s/%s'                       % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sOutDir, exist_ok=True)

    sScript         = 'date; %s mt_correction'          % sTempScript
    sScript        += ' %s %s %s; date'                 % (sInDir, sOutDir,  nJobs)
    sScript         = 'echo "%s" | qsub -cwd -j y -o %s -q %s -N %s' \
                      % (sScript, sLogDir, sQueue, sJobName)
    if bTestRun: print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system(sScript)
    #if END:
#def END: qsub_multipletest


def mt_correction (sInDir, sOutDir, nJobs):
    print_start_msg('Multiple test correction')
    list_cHotSpots = []
    for nJobIndex in range(int(nJobs)):
        ## MD
        sInFile         = '%s/%04d.dat' % (sInDir, nJobIndex+1)
        InFile          = open(sInFile, 'rb')
        list_cHotSpots  += pickle.load(InFile)
        InFile.close()
    #loop END: nJobIndex
    list_fPvals  = [[cHotSpot.fPvalue, cHotSpot.nClusterID] for cHotSpot in list_cHotSpots]

    list_fQvals  = get_qval_list(list_fPvals)
    dict_fQvals  = dict(list_fQvals)

    list_sOutput = []
    for cHotSpot in list_cHotSpots:
        cHotSpot.fQvalue = dict_fQvals[cHotSpot.nClusterID]
        list_sOutput.append(cHotSpot)
    #loop END: cHotSpot

    OutFile      = open('%s/filtered_hotspots.dat' % (sOutDir), 'wb')
    pickle.dump(list_sOutput, OutFile)
    OutFile.close()

    list_sOutput = sorted(list_sOutput, key=lambda c:(c.sChrID.replace('chr',''), c.nStartPos))
    OutFile      = open('%s/filtered_hotspots_stats.txt' % (sOutDir), 'w')
    for cHotspot in list_sOutput:
        sPosList  = ','.join([str(nPos) for nPos in cHotspot.list_nPos])
        sOut      = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                   (cHotspot.sChrID, cHotspot.nStartPos, cHotspot.nEndPos, cHotspot.nClusterID,
                    cHotspot.nSize, cHotspot.nVarCnt, len(cHotspot.list_sPatIDs), cHotspot.nPosCnt,
                    cHotspot.fQvalue, cHotspot.fBackProb, cHotspot.sLocale, cHotspot.sGeneSym,
                    cHotspot.sGeneID, sPosList)
        OutFile.write(sOut)
    #loop END: cHotspot
    OutFile.close()
    print_done_msg('Multiple test correction')
#def END: mt_correction


def basic_stats_hotspots (nClustDist, nBackWinSize):
    print(green('STEP 5: Survey statistically significant (P<0.05) hotspots', 1))
    sHotspotDir     = '%s/04_mtcorrected/%sbps/%s'     % (sBASE_DIR, nClustDist, nBackWinSize)
    sHotspotFile    = '%s/filtered_hotspots.dat'       % sHotspotDir

    InFile          = open(sHotspotFile, 'rb')
    list_cHotspots  = pickle.load(InFile)
    InFile.close()

    list_cHotspots  = sorted(list_cHotspots, key=lambda c:(c.sChrID.replace('chr',''), c.nStartPos))

    #V-S Check
    if not list_cHotspots:
        sys.exit('Invalid list : list_cHotspots size= %d' % len(list_cHotspots))

    dict_sLocale  = {}
    for cHotspot in list_cHotspots:
        sKey = dict_sGENIC_alt[cHotspot.sLocale]
        if sKey not in dict_sLocale:
            dict_sLocale[sKey] = 0
        dict_sLocale[sKey] += 1
    #loop END: cHotspot

    for sLocale in dict_sLocale:
        print(sLocale, dict_sLocale[sLocale])
#def END: basic_stats_hotspots


def qsub_get_rbp_variants (list_cClipData, nClustDist, nBackWinSize, sTempScript,  sQueue, bTestRun):
    print(green('STEP 6: Overlap with Bedtools', 1))
    sJobName      = 'Jinman.Intersect.RBP.%s'        % nBackWinSize
    sInDir        = '%s/04_mtcorrected/%sbps/%s'     % (sBASE_DIR, nClustDist, nBackWinSize)
    sClipDir      = '%s/sorted_bed'                  % sECLIP_DIR
    sLogDir       = '%s/log/%s/%s'                   % (sBASE_DIR, sJobName, sTIME_STAMP)
    sOutDir       = '%s/05_rbp_intersect/%sbps/%s'   % (sBASE_DIR, nClustDist, nBackWinSize)
    os.makedirs(sOutDir, exist_ok=True)

    sHotspotFile  = '%s/filtered_hotspots_stats.txt' % (sInDir)

    for cClipData in list_cClipData:

        sScript = '%s get_rbp_variants %s %s %s %s' \
                    % (sTempScript, sHotspotFile, sClipDir, sOutDir, cClipData.sAccID)

        sScript = 'echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s' \
                    % (sScript, sLogDir, sQueue, sJobName, cClipData.sAccID)

        if bTestRun: print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system(sScript)
    #loop END: cClipData
#def END: qsub_get_rbp_variants


def qsub_get_rbp_variants_pt2 (list_cClipData, nClustDist, nBackWinSize, nBufferSize, sFileTag, sTempScript,  sQueue, bTestRun):
    print(green('STEP 9-2: Overlap with Bedtools Part 2', 1))
    sJobName      = 'Jinman.Intersect.RBP.Pt2.%s'              % nBackWinSize
    sInDir        = '%s/08_rbps_wInters/%sbps/%s'              % (sBASE_DIR, nClustDist, nBackWinSize)
    sClipDir      = '%s/01_basicstats/ClipData_buffered/%sbps' % (sBASE_DIR, nBufferSize)
    sLogDir       = '%s/log/%s/%s'                             % (sBASE_DIR, sJobName, sTIME_STAMP)
    sOutDir       = '%s/08_rbps_wInters/%sbps/%s/out%s'        % (sBASE_DIR, nClustDist, nBackWinSize, nBufferSize)
    os.makedirs(sOutDir, exist_ok=True)

    sHotspotFile  = '%s/rbphotspot_forBuffer_%s.txt' % (sInDir, sFileTag)

    for cClipData in list_cClipData:

        sScript = '%s get_rbp_variants %s %s %s %s' \
                    % (sTempScript, sHotspotFile, sClipDir, sOutDir, cClipData.sAccID)

        sScript = 'echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s' \
                    % (sScript, sLogDir, sQueue, sJobName, cClipData.sAccID)

        if bTestRun: print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system(sScript)
    #loop END: cClipData
#def END: qsub_get_rbp_variants_pt2


def get_rbp_variants (sHotspotFile, sClipDir, sOutDir, sAccID):
    print_start_msg('Intersecting BedFile', sBASE_DIR)
    sBedFile     = '%s/%s.bed'           % (sClipDir, sAccID)
    if not os.path.isfile(sBedFile): sys.exit('File Not Found %s' % sBedFile)

    sScript      = 'bedtools intersect -a %s -b %s -wao'  % (sHotspotFile, sBedFile)
    sStdOut      = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout
    list_sOutput = []

    for sReadLine in sStdOut:
        #Standard out format
        '''
        0   chr1                        sChrID                  Hotspot Info
        1   1334663                     nHotspotStart
        2   1334663                     nHotSpotEnd
        3   1497                        nHotspotID
        4   1                           nSize
        5   9                           nVarCnt
        6   9                           nPatCnt
        7   0                           nPosCnt
        8   8.178838834508808e-24       fQvalue
        9   2.05-06                     fBackProb
        10  exonic                      sLocale
        11  CCNL2                       sGeneSym
        12  ENSG00000221978             sGeneID
        13  chr1                        sChrID                  RBP Info
        14  1334580                     nRBPStart
        15  1334663                     nRBPEnd
        16  GRSF1_HepG2_rep02           sRBPID
        17  200                         nSomeNumber
        18  -                           sStrand
        19  -0.754152666559453          fScore
        20  0                           NULL
        21  -1                          NULL
        22  -1                          NULL
        23  1                           sOverlap
        '''
        sReadLine         = str(sReadLine, 'UTF-8')
        list_sColumn      = sReadLine.strip('\n').split('\t')

        cRBP              = cRBP_Cluster()
        cRBP.sChrID       = list_sColumn[0]
        cRBP.nStartPos    = int(list_sColumn[1])
        cRBP.nEndPos      = int(list_sColumn[2])
        cRBP.nClusterID   = list_sColumn[3]
        cRBP.nSize        = int(list_sColumn[4])
        cRBP.nVarCnt      = int(list_sColumn[5])
        cRBP.nPatCnt      = int(list_sColumn[6])
        cRBP.nPosCnt      = int(list_sColumn[7])
        cRBP.fQvalue      = float(list_sColumn[8])
        cRBP.fBackProb    = float(list_sColumn[9])
        cRBP.sLocale      = list_sColumn[10]
        cRBP.sGeneSym     = list_sColumn[11]
        cRBP.sGeneID      = list_sColumn[12]
        cRBP.list_nPos    = [int(sPos) for sPos in list_sColumn[13].split(',')]
        cRBP.nRBPStart    = int(list_sColumn[15])
        cRBP.nRBPEnd      = int(list_sColumn[16])
        cRBP.sRBPID       = list_sColumn[17]
        cRBP.nOverlap     = int(list_sColumn[24]) if cRBP.nSize >= 1  else 0
        if cRBP.sRBPID == '.': cRBP.nOverlap   = 0
        '''
        if cRBP.sRBPID != '.':
            sOut = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
               % (cRBP.nClusterID, cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos,
                  cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nPosCnt, cRBP.fQvalue, cRBP.fQvalue,
                  dict_sGENIC[cRBP.sLocale], cRBP.sGeneSym, cRBP.sGeneID, cRBP.list_nPos,
                  cRBP.nRBPStart, cRBP.nRBPEnd, cRBP.sRBPID.split('_')[0], cRBP.nOverlap)
            print(sOut[:-1])
            sys.exit()
        '''
        list_sOutput.append(cRBP)
    #loop END: sReadLine

    sOutFile     = '%s/%s.intersect.dat' % (sOutDir, sAccID)
    OutFile      = open(sOutFile, 'wb')
    pickle.dump(list_sOutput, OutFile)
    OutFile.close()
    print_done_msg('Annotating BedFile', sBASE_DIR)
#def END: get_rbp_variants


def basic_stats_clusters_rbps (list_cClipData, nClustDist, nBackWinSize, fRecurrency):
    print(green('STEP 7-1: Basic Statistics of Clusters-RBPs P<0.05', 1))
    sOutDir = '%s/06_rbp_hotspots/%sbps/%s' % (sBASE_DIR, nClustDist, nBackWinSize)
    os.makedirs(sOutDir, exist_ok=True)

    ## Load COSMIC data ##
    #dict_cCosGene = load_cosmic()

    ## Load CPDB data ##
    #dict_sCPDBGene = load_CPDB()

    list_sLocale  = ['exonic', '3UTR', '5UTR', 'intergenic', 'intronic', 'splicing', 'upstream', 'downstream']
    sInDir        = '%s/05_rbp_intersect/%sbps/%s' % (sBASE_DIR, nClustDist, nBackWinSize)
    list_sAllRBPs = []
    for cClipData in list_cClipData:
        sInFile         = '%s/%s.intersect.dat' % (sInDir, cClipData.sAccID)
        #list_sIntersect = parse_intersect_file(sInFile)
        InFile          = open(sInFile, 'rb')
        list_sIntersect = pickle.load(InFile)
        InFile.close()
        list_sAllRBPs  += list_sIntersect
    #loop END: cClipData

    print('list_sAllRBPs', len(list_sAllRBPs))

    dict_nClusterID  = {}
    for cRBP in list_sAllRBPs:

        #if cRBP.fQvalue > 0.05: continue

        #if cRBP.sGeneSym not in dict_sCPDBGene: #CPDB ovelap
            #continue

        #if cRBP.sLocale == 'exonic': continue #Filter exonic

        #if fRecurrency > 0: #Recurrency Filter ON
            #if cRBP.nPatCnt <= int(nSAMPLE_CNT * 0.10): continue
        #if END:
        sKey = cRBP.nClusterID

        if sKey not in dict_nClusterID:
            dict_nClusterID[sKey] = []
        dict_nClusterID[sKey].append(cRBP)
    #loop END: cRBP
    print('Number of RBP-Cluster Interactions', len(dict_nClusterID))

    list_sFinalRBPs = []
    list_nClusterID = list(dict_nClusterID.keys())
    for nClusterID in list_nClusterID:

        dict_sRBPID = {}

        for cRBP in dict_nClusterID[nClusterID]:

            sKey = cRBP.sRBPID.split('_')[0]
            if sKey not in dict_sRBPID:
                dict_sRBPID[sKey] = []
            dict_sRBPID[sKey].append(cRBP)
        #loop END: cRBP

        for sRBPID in dict_sRBPID:

            if len(dict_sRBPID[sRBPID]) == 1: continue
            list_sFinalRBPs.append(dict_sRBPID[sRBPID][0])
            #loop END: cRBP
        #loop END: sRBPID

    #loop END: nClusterID
    print('FINAL Number of RBP-Cluster Interactions', len(list_sFinalRBPs))
    sys.exit()

    #OutFile = open('%s/final_rbphotspot.txt' % sOutDir, 'w')
    dict_sLocale    = {}
    for i, cRBP in enumerate(list_sFinalRBPs):

        sKey = dict_sGENIC[cRBP.sLocale]

        if sKey not in dict_sLocale:
            dict_sLocale[sKey] = []
        dict_sLocale[sKey].append(cRBP)

        sOut = ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                (i, cRBP.nClusterID, cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos,
                cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nPosCnt, cRBP.fQvalue,
                dict_sGENIC[cRBP.sLocale], cRBP.sGeneSym, cRBP.sRBPID.split('_')[0], cRBP.nOverlap))
        print(sOut[:-1])
        #OutFile.write(sOut)
    #loop END: cRBP
    #OutFile.close()

    for sLocale in dict_sLocale:
        print(sLocale, len(dict_sLocale[sLocale]))
    #loop END: sLocale
#def END: basic_stats_clusters_rbps


def compile_rbp_data (list_cClipData, nClustDist, nBackWinSize, sFileTag):
    print(green('STEP 7: Compile RBP data', 1))
    sOutDir       = '%s/06_rbp_hotspots/%sbps/%s'  % (sBASE_DIR, nClustDist, nBackWinSize)
    os.makedirs(sOutDir, exist_ok=True)

    sInDir        = '%s/05_rbp_intersect/%sbps/%s' % (sBASE_DIR, nClustDist, nBackWinSize)
    list_sAllRBPs = []
    for cClipData in list_cClipData:
        sInFile         = '%s/%s.intersect.dat'    % (sInDir, cClipData.sAccID)
        #list_sIntersect = parse_intersect_file(sInFile)
        InFile          = open(sInFile, 'rb')
        list_sIntersect = pickle.load(InFile)
        InFile.close()
        list_sAllRBPs  += list_sIntersect
    #loop END: cClipData
    print('Loaded RBP data', len(list_sAllRBPs))

    list_cRBP = consolidate_rbplist(list_sAllRBPs)
    print('FINAL Number of RBP-Cluster All Interactions', len(list_cRBP))

    # Pickle Dumpage
    sOutFile = '%s/rbphotspot_%s.dat' % (sOutDir, sFileTag)
    output_rbp_data(sOutFile, list_cRBP)
#def END: compile_rbp_data


def load_cosmic ():
    sCosmicFile    = '%s/cosmic.v74.txt'         % sREF_DIR
    dict_sCosGene  = {}
    list_cCosmic   = cCos_parse_cosmic_consensus (sCosmicFile)

    for cCos in list_cCosmic:
        if cCos.sGeneName not in dict_sCosGene:
            dict_sCosGene[cCos.sGeneName] = []
        dict_sCosGene[cCos.sGeneName].append(cCos)
    #loop END: cCos
    print('Cosmic Genes', len(dict_sCosGene))

    return dict_sCosGene
#def END: load_cosmic


def load_CPDB ():
    ## Load CPDB data ##
    sCPDBFile       = '%s/CPDB_pathways.txt'        % sREF_DIR
    dict_sCPDBGene  = parse_pathways (sCPDBFile, 0)

    return dict_sCPDBGene
#def END: load_CPDB


def parse_intersect_file (sInFile):

    #VS-Check
    if not os.path.isfile(sInFile):
        sys.exit('File Not Found %s' % sInFile)

    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:
        list_sColumn    = sReadLine.strip('\n').split('\t')

        cRBP            = cRBP_Cluster()
        cRBP.sChrID     = list_sColumn[0]
        cRBP.nStartPos  = int(list_sColumn[1])
        cRBP.nEndPos    = int(list_sColumn[2])
        cRBP.nClusterID = int(list_sColumn[3])
        cRBP.nSize      = int(list_sColumn[4])
        cRBP.nVarCnt    = int(list_sColumn[5])
        cRBP.nPatCnt    = int(list_sColumn[6])
        cRBP.nCovered   = int(list_sColumn[7])
        cRBP.fQvalue    = float(list_sColumn[8])
        cRBP.fBackProb  = float(list_sColumn[9])
        cRBP.sLocale    = list_sColumn[10]
        cRBP.sGeneSym   = list_sColumn[11]
        cRBP.sGeneID    = list_sColumn[12]
        cRBP.sRBPID     = list_sColumn[13]
        cRBP.nOverlap   = int(list_sColumn[14])
        cRBP.nTotalSize = int(list_sColumn[15])
        if cRBP.sRBPID == '.': cRBP.nOverlap = 0

        list_sOutput.append(cRBP)
    #loop END: sReadLine
    InFile.close()
    #VS-Check
    if not list_sOutput:
        sys.exit('Invalid List : parse_intersect_file : list_sOutput size= %d' % len(list_sOutput))

    #for e in dict_sTest:
        #print(e, dict_sTest[e])
    #sys.exit()
    return list_sOutput
#def END: parse_intersect_file


def parse_pathways (sInFile, bCancer):
    print_start_msg('Parsing CPDB File', sInFile)
    dict_sGeneSym   = {}
    InFile          = open(sInFile, 'r')
    dict_sPathways  = {}
    for sReadLine in InFile:
        if sReadLine.startswith('pathway'): continue
        list_sColumn = sReadLine.strip('\n').split('\t')

        sPathName       = list_sColumn[0].lower()
        sExternalID     = list_sColumn[1]
        sSource         = list_sColumn[2]
        list_sGeneSym   = [sGene.upper() for sGene in  list_sColumn[3].split(',')]

        if bCancer:
            if 'cancer' not in sPathName: continue ## Filter out non-cancer related pathways

        if sPathName not in dict_sPathways:
            dict_sPathways[sPathName] = []
        dict_sPathways[sPathName] = list_sGeneSym

        for sGeneSym in list_sGeneSym:
            sKey = sGeneSym.upper()
            if sKey not in dict_sGeneSym:
                dict_sGeneSym[sKey] = []
            dict_sGeneSym[sKey].append(sPathName)
        #loop END: sGeneSym
    #loop END: sReadLine
    InFile.close()

    #for sPathName in dict_sPathways:
        #print(sPathName)

    #sys.exit()

    #V-S Check:
    if not dict_sGeneSym:
        sys.exit('Invalid Dictionary : parse_pathways : dict_sGeneSym size= %d' % len(dict_sGeneSym))

    if bCancer:
        print('CPDB Cancer Genes', len(dict_sGeneSym))

    return dict_sGeneSym
#def END: parse_pathways


def output_rbp_data (sOutFile, list_cOutput):
    OutFile = open(sOutFile, 'wb')
    pickle.dump(list_cOutput, OutFile)
    OutFile.close()
    OutFile = open(sOutFile.replace('.dat', '.txt'), 'w')
    for i, cRBP in enumerate(list_cOutput):
        sOut = ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                i, cRBP.nClusterID, cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos,
                cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nPosCnt, cRBP.fQvalue,
                dict_sGENIC[cRBP.sLocale], cRBP.sGeneSym, cRBP.sRBPID.split('_')[0], cRBP.nOverlap))
        OutFile.write(sOut)
    #loop END: cRBP
    OutFile.close()
#def END: output_rbp_data


def output_rbp_data_v2 (list_cOutput):

    sHeader         = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
                      % ('ChrID','Start','End','Size','Patients','RNAseq','fQvalue','Region', 'TargetGene',
                         'RBPs', 'OverlapSize', 'wHot', 'woHot')
    print(sHeader)

    for cRBP in list_cOutput:
        sOut              = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                            (cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos, cRBP.nSize, cRBP.nPatCnt,
                             cRBP.nPosCnt, cRBP.fQvalue, cRBP.sLocale, cRBP.sGeneSym,
                             len(cRBP.list_sRBPID), cRBP.nOverlap, cRBP.fExprDeviation_M,
                             cRBP.fExprDeviation_NM)
        print(sOut[:-1])
    #loop END: cRBP
#def END: output_rbp_data


def qsub_preprocess_fpkm (sTempScript, sQueue, bTestRun):
    print(green('STEP 8-1: Preprocess TCGA FPKM files', 1))
    # Load TCGA FPKM File List
    sFPKMFile       = '%s/%s/%s/%s_rnaseq-UQ.txt' % (sBASE_DIR, sSTUDY, sSEQDATA, sSTUDY)
    dict_sFPKMFiles = load_fpkm_list (sFPKMFile)

    sInDir   = '%s/%s/FPKMfiles/UQ_norm'          % (sBASE_DIR, sSTUDY)
    sOutDir  = '%s/%s/FPKMfiles/UQ_norm_recode'   % (sBASE_DIR, sSTUDY)
    sJobName = 'Jinman.Preprocess.FPKM'
    sLogDir  = '%s/log/%s/%s'                     % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    for sPatID in dict_sFPKMFiles:
        sInFile_N  = '%s/%s'                     % (sInDir, dict_sFPKMFiles[sPatID]['N'].replace('.gz', ''))
        sInFile_T  = '%s/%s'                     % (sInDir, dict_sFPKMFiles[sPatID]['T'].replace('.gz', ''))
        sOutFile_N = '%s/%s.%s.txt'              % (sOutDir, sPatID, 'N')
        sOutFile_T = '%s/%s.%s.txt'              % (sOutDir, sPatID, 'T')
        sScript    = 'date; %s preprocess_fpkm ' % sTempScript
        sScript   += '%s %s %s %s %s; date'      % (sInFile_N, sInFile_T , sOutFile_N, sOutFile_T, fTOP_EXPR)
        sScript    = 'echo "%s" | qsub -cwd -j y -o %s.log -q %s -N %s' \
                      % (sScript, sPatID, sQueue, sJobName)
        if bTestRun: print(sScript)
        else:
            os.system(sScript)
        #if END:
    #loop END: sPatID
#def END: preprocess_fpkm


def preprocess_fpkm (sInFile_N, sInFile_T, sOutFile_N, sOutFile_T, fTOP_EXPR):
    # Load Gencode Reference Data
    sGencodeFile    = '%s/gencode.v22.annotation.gtf' % (sREF_DIR)
    dict_sGencode   = load_reference_gencode (sGencodeFile)

    dict_sOutput_N  = load_raw_fpkm_file(sInFile_N, dict_sGencode)
    dict_sOutput_T  = load_raw_fpkm_file(sInFile_T, dict_sGencode)
    list_sTopExpr   = [[sGeneSym] + dict_sOutput_N[sGeneSym] for sGeneSym in dict_sOutput_N]
    list_sTopExpr   = sorted(list_sTopExpr, key=lambda e: e[2], reverse=True)[:int(len(dict_sOutput_N)* float(fTOP_EXPR))]
    write_processed_fpkm_file(sOutFile_N, dict_sOutput_N, list_sTopExpr)
    write_processed_fpkm_file(sOutFile_T, dict_sOutput_T, list_sTopExpr)
    #loop END: i, sPatID
#def END: preprocess_fpkm


def load_raw_fpkm_file (sInFile, dict_sGencode):
    InFile        = open(sInFile, 'r')
    dict_sOutput  = {}
    for sReadLine in InFile:
        list_sColumn = sReadLine.strip('\n').split('\t')
        sGeneID      = list_sColumn[0]
        fFPKM        = float(list_sColumn[1])

        #if fFPKM <= 0: continue # No expression
        try:             sGeneSym = dict_sGencode[sGeneID]
        except KeyError: continue # GeneID not found in Gencode

        if sGeneSym not in dict_sOutput:
            dict_sOutput[sGeneSym] = ''
        dict_sOutput[sGeneSym] = [sGeneID.split('.')[0], fFPKM]

    #loop END: sReadLine
    InFile.close()

    #VS-Check
    if not dict_sOutput:
        sys.exit('Invalid List : load_raw_fpkm_file : dict_sOutput size= %d' % len(dict_sOutput))

    return dict_sOutput
#def END: load_raw_fpkm_file


def write_processed_fpkm_file (sOutFile, dict_sOutput, list_sTopExpr):
    OutFile       = open(sOutFile, 'w')
    for sGeneSym, sGeneID, fFPKM_N in list_sTopExpr:
        sGeneID, fFPKM = dict_sOutput[sGeneSym]
        sOut           = '%s\t%s\t%s\n' % (sGeneSym, sGeneID.split('.')[0], fFPKM)
        OutFile.write(sOut)
    #loop END: sGeneID, fFPKM
    OutFile.close()
#def END: write_processed_fpkm_file


def evalute_expression_data (nClustDist, nBackWinSize):
    print(green('STEP 8: Expression Data Comparison', 1))
    # Open original hotspot file (for patient id list)
    sHotspotDir   = '%s/04_mtcorrected/%sbps/%s'        % (sBASE_DIR, nClustDist, nBackWinSize)
    sHotspotFile  = '%s/filtered_hotspots.%s.dat'       % (sHotspotDir, 0)
    InFile        = open(sHotspotFile, 'rb')
    list_cHotspot = pickle.load(InFile)
    InFile.close()

    # Open final rbp-hotspot file
    sInDir        = '%s/06_rbp_hotspots/%sbps/%s' % (sBASE_DIR, nClustDist, nBackWinSize)
    sInFile       = '%s/rbphotspot.dat'           % sInDir
    InFile        = open(sInFile, 'rb')
    list_cRBP     = pickle.load(InFile)
    InFile.close()

    list_cRBP     = get_final_list (list_cHotspot, list_cRBP)

    sFPKMDir      = '%s/%s/FPKMfiles/UQ_norm_recode'    % (sBASE_DIR, sSTUDY)

    for cRBP in list_cRBP:
        list_fFPKM_T = []
        list_fFPKM_N = []
        list_fFC     = []
        for sPatID in cRBP.list_sCovered:
            sFPKMFile_T = '%s/%s.T.txt' % (sFPKMDir, sPatID)
            sFPKMFile_N = '%s/%s.N.txt' % (sFPKMDir, sPatID)

            fFPKM_T     = load_fpkm_v2(sFPKMFile_T, cRBP.sGeneSym)
            fFPKM_N     = load_fpkm_v2(sFPKMFile_N, cRBP.sGeneSym)

            list_fFPKM_T.append(fFPKM_T)
            list_fFPKM_N.append(fFPKM_N)
            list_fFC.append(math.log(fFPKM_T / fFPKM_N, 2))
        #loop END: sPatID

        sOut = ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'
                % (cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos, cRBP.nSize,
                   cRBP.nPatCnt, cRBP.nPosCnt, cRBP.fQvalue, cRBP.sLocale,
                   cRBP.sGeneSym, len(cRBP.list_sRBPID), cRBP.nOverlap, np.mean(list_fFC), stats.sem(list_fFC)))
        print(sOut[:-1])
    #loop END: cRBP
#def END: evalute_expression_data


def evalute_expression_data_v2 (nClustDist, nBackWinSize, fRecurrency):
    print(green('STEP 8: Expression Data Analysis', 1))

    nBins           = 4
    # Load TCGA FPKM File List
    sFPKMDir        = '%s/%s/FPKMfiles/UQ_norm_recode'    % (sBASE_DIR, sSTUDY)
    sFPKMFile       = '%s/%s/%s_rnaseq-UQ.txt'            % (sBASE_DIR, sSTUDY, sSTUDY)
    dict_sFPKMFiles = load_fpkm_list (sFPKMFile)

    # Open original hotspot file (for patient id list)
    sHotspotDir     = '%s/04_mtcorrected/%sbps/%s'        % (sBASE_DIR, nClustDist, nBackWinSize)
    sHotspotFile    = '%s/filtered_hotspots.%s.dat'       % (sHotspotDir, 0)
    InFile          = open(sHotspotFile, 'rb')
    list_cHotspot   = pickle.load(InFile)
    InFile.close()

    # Open final rbp-hotspot file
    sInDir          = '%s/06_rbp_hotspots/%sbps/%s' % (sBASE_DIR, nClustDist, nBackWinSize)
    sInFile         = '%s/rbphotspot_recurrent.dat' % sInDir
    InFile          = open(sInFile, 'rb')
    list_cRBP       = pickle.load(InFile)
    InFile.close()

    # Final non-redundant RBP list
    list_cRBP       = get_final_list (list_cHotspot, list_cRBP)
    dict_sExpr      = assign_fpkm_values (sFPKMDir, list_cRBP, dict_sFPKMFiles)

    list_sOutput  = []
    for cRBP in list_cRBP:
        #if cRBP.sGeneSym not in set(list(dict_sCPDBGene.keys())): continue
        if fRecurrency > 0: #Recurrency Filter ON
            if cRBP.nPatCnt <= int(nSAMPLE_CNT * 0.1): continue

        # Mutated Samples
        list_fFC_M  = [dict_sExpr[sPatID][cRBP.sGeneSym] for sPatID in cRBP.list_sCovered]
        if None in list_fFC_M: continue

        # Non-mutated Samples
        list_fFC_NM = [dict_sExpr[sPatID][cRBP.sGeneSym] for sPatID in dict_sExpr if sPatID not in set(cRBP.list_sCovered)]
        if None in list_fFC_NM: continue

        cRBP.fExprDeviation_M  = calculate_expr_deviation (list_fFC_M, list_fFC_NM, 'M')
        cRBP.fExprDeviation_NM = calculate_expr_deviation (list_fFC_M, list_fFC_NM, 'NM')

        list_sOutput.append(cRBP)
    #loop END: cRBP
    list_sOutput       = sorted(list_sOutput, key= lambda c:c.fExprDeviation_M - c.fExprDeviation_NM, reverse=True)

    output_rbp_data_v2(list_sOutput)
    sys.exit()
    print('Final RBP list', len(list_sOutput))

    list_sBinned    = [[] for i in range(nBins)]

    assign_bin(list_sOutput, list_sBinned)

    for i, list_cRBP_binned in enumerate(list_sBinned):
        list_fExpr_M  = []
        list_fExpr_NM = []

        for cRBP in list_cRBP_binned:
            list_fExpr_M.append(cRBP.fExprDeviation_M)
            list_fExpr_NM.append(cRBP.fExprDeviation_NM)
        #loop END: cRBP
        Z, fPvalue = stats.ranksums(list_fExpr_M, list_fExpr_NM)
        print(i, Z, fPvalue)
    #loop END: list_cRBP_binned
#def END: evalute_expression_data


def assign_bin(list_sOutput, list_sBinned):

    nTotalOutput    = len(list_sOutput)
    nBins           = len(list_sBinned)

    for i, cRBP in enumerate(list_sOutput):
        nBin = int(i/nTotalOutput * nBins)
        list_sBinned[nBin].append(cRBP)
    return list_sBinned
#def END: assign_bin


def assign_fpkm_values (sFPKMDir, list_cRBP, dict_sFPKMFiles):
    dict_sTargetGenes = {}
    for cRBP in list_cRBP:
        if cRBP.sGeneSym not in dict_sTargetGenes:
            dict_sTargetGenes[cRBP.sGeneSym] = []
        dict_sTargetGenes[cRBP.sGeneSym].append(cRBP)
    #loop END: sGenSym

    #list_sGrepGenes   = ['\\b%s\\b' % sGeneSym for sGeneSym in dict_sTargetGenes]
    dict_sOutput      = {}

    for sPatID in dict_sFPKMFiles:

        if sPatID not in dict_sOutput:
            dict_sOutput[sPatID] = ''
        dict_sOutput[sPatID] = {}

        sFPKMFile   = '%s/%s.T.txt' % (sFPKMDir, sPatID)
        dict_fFPKM  = load_fpkm(sFPKMFile)

        for sGeneSym in dict_sTargetGenes:
            try:             fFPKM = dict_fFPKM[sGeneSym]
            except KeyError: fFPKM = None

            if sGeneSym not in dict_sOutput[sPatID]:
                dict_sOutput[sPatID][sGeneSym] = 0.0
            dict_sOutput[sPatID][sGeneSym] = fFPKM
        #loop END: sGeneSym
    #loop END: sPatID
    #VS-Check
    if not dict_sOutput:
        sys.exit('Invalid Dict : assign_fpkm_values : dict_sOutput size= %d' % len(dict_sOutput))

    return dict_sOutput
#def END: assign_fpkm_values


def load_fpkm (sInFile):

    dict_sOutput = {}
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:
        list_sColumn = sReadLine.strip('\n').split('\t')
        sGeneSym     = list_sColumn[0]
        sGeneID      = list_sColumn[1]
        fFPKM        = float(list_sColumn[2])

        if sGeneSym not in dict_sOutput:
            dict_sOutput[sGeneSym] = []
        dict_sOutput[sGeneSym].append(fFPKM)
    #loop END:
    InFile.close()

    dict_fFPKM = {}
    for sGeneSym in dict_sOutput:

        fFPKM = np.mean(list(filter(None, dict_sOutput[sGeneSym]))) if list(filter(None, dict_sOutput[sGeneSym])) else 0
        if sGeneSym not in dict_fFPKM:
            dict_fFPKM[sGeneSym] = ''
        dict_fFPKM[sGeneSym] = fFPKM if fFPKM > 0 else 0.0
    #loop END: sGeneSym

    #VS-Check
    if not dict_fFPKM:
        sys.exit('Invalid Dict : load_fpkm : dict_fFPKM size= %d' % len(dict_fFPKM))

    return dict_fFPKM
#def END: load_fpkm


def load_fpkm_v2 (sInFile, sTargetGeneID):

    sScript = 'egrep \"\\b%s\\b\" %s' % (sTargetGeneID, sInFile)
    sStdOut = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout
    fFPKM   = '-1'
    for sReadLine in sStdOut:
        sGeneSym, sGeneID, fFPKM = (str(sReadLine, 'UTF-8').strip('\n').split('\t'))
    #loop END: sReadLine
    return float(fFPKM)
#def END: load_fpkm_v2


def get_final_list (list_cHotspot, list_cRBP):
    print('list_cHotspot', len(list_cHotspot))
    dict_sCovered = {}
    for cHotspot in list_cHotspot:
        sKey      = '%s.%s.%s' % (cHotspot.nClusterID, cHotspot.sChrID, len(cHotspot.list_sCovered))

        if sKey not in dict_sCovered:
            dict_sCovered[sKey] = ''
        dict_sCovered[sKey] = cHotspot
    #loop END: sKey

    print('list_cRBP', len(list_cRBP))
    dict_nClusterID = {}
    for cRBP in list_cRBP:

        if cRBP.nCovered == 0: continue  # Filter out interactions with no transcriptome data

        if cRBP.nClusterID not in dict_nClusterID:
            dict_nClusterID[cRBP.nClusterID] = []
        dict_nClusterID[cRBP.nClusterID].append(cRBP)
    #loop END: cRBP


    list_cOutput = []
    for nClusterID in dict_nClusterID:

        cOutput               = cRBP_Cluster()
        cOutput.nClusterID    = nClusterID
        cOutput.sChrID        = [cRBP.sChrID for cRBP in dict_nClusterID[nClusterID]][0]
        cOutput.nStartPos     = [cRBP.nStartPos for cRBP in dict_nClusterID[nClusterID]][0]
        cOutput.nEndPos       = [cRBP.nEndPos for cRBP in dict_nClusterID[nClusterID]][0]
        cOutput.nSize         = [cRBP.nSize for cRBP in dict_nClusterID[nClusterID]][0]
        cOutput.nVarCnt       = [cRBP.nVarCnt for cRBP in dict_nClusterID[nClusterID]][0]
        cOutput.nPatCnt       = [cRBP.nPatCnt for cRBP in dict_nClusterID[nClusterID]][0]
        cOutput.nCovered      = [cRBP.nCovered for cRBP in dict_nClusterID[nClusterID]][0]
        cOutput.fQvalue       = [cRBP.fQvalue for cRBP in dict_nClusterID[nClusterID]][0]
        cOutput.sLocale       = [dict_sGENIC[cRBP.sLocale] for cRBP in dict_nClusterID[nClusterID]][0]
        cOutput.sGeneSym      = [cRBP.sGeneSym for cRBP in dict_nClusterID[nClusterID]][0]
        cOutput.list_sRBPID   = [cRBP.sRBPID.split('_')[0] for cRBP in dict_nClusterID[nClusterID]]
        cOutput.nOverlap      = [cRBP.nOverlap for cRBP in dict_nClusterID[nClusterID]][0]

        sKey                  = '%s.%s.%s' % (nClusterID, cOutput.sChrID, cOutput.nCovered)
        cOutput.list_sCovered = dict_sCovered[sKey].list_sCovered
        cOutput.sGeneID       = dict_sCovered[sKey].sGeneID

        list_cOutput.append(cOutput)
    #loop END: nClusterID

    #VS-Check
    if not list_cOutput:
        sys.exit('Invalid List : get_final_list : list_cOutput size= %d' % len(list_cOutput))
    print('Final RBP list', len(list_cOutput))
    return list_cOutput
#def END: get_final_list


def calculate_expr_deviation (list_fFC_M, list_fFC_NM, bOutput):
    fMeanExpr_all  = np.mean(list_fFC_M + list_fFC_NM)
    list_fFC       = list_fFC_M if bOutput == 'M' else list_fFC_NM
    sGroupSize     = len(list_fFC)
    fExprDeviation = sum([abs(fFC - fMeanExpr_all) / fMeanExpr_all for fFC in list_fFC]) / sGroupSize
    return fExprDeviation
#def END: calculate_expr_deviation


def determine_miRNA_targetsites (sGenome, nClustDist, nBackWinSize, nBufferSize, sFileTag):
    print(green('STEP 8: Compare miRNA target sites', 1))

    '''Load 3UTRs n=18625'''
    bCapture        = 0
    dict_c3UTRs     = c3UTRs_load_3UTRS(s3UTR_FILE[sGenome])
    list_sCapture   = parse_capture_bedfile(sCAPTUREFILE)
    dict_cRefSeq    = cRef_parse_refflat_line(sREFSEQFILE, list_sCapture, bCapture)

    '''Assign 3UTR genomic position info'''
    annotate_3UTRs (dict_c3UTRs, dict_cRefSeq)

    sOutDir         = '%s/07_rbps_wMirTargets/%sbps/%s'   % (sBASE_DIR, nClustDist, nBackWinSize)
    os.makedirs(sOutDir, exist_ok=True)

    # Load TCGA FPKM File List
    sFPKMDir        = '%s/%s/FPKMfiles/UQ_norm_recode'    % (sBASE_DIR, sSTUDY)
    sFPKMFile       = '%s/%s/%s_rnaseq-UQ.txt'            % (sBASE_DIR, sSTUDY, sSTUDY)
    dict_sFPKMFiles = load_fpkm_list (sFPKMFile)
    dict_sExpr      = get_expr_data (sFPKMDir, dict_sFPKMFiles)

    # Obtain top miRNAs target sites in HepG2
    list_cMir       = load_top_mirnas (dict_sMIRFILE[sCELLS])
    dict_sTargets   = determine_site_type_cnt (dict_c3UTRs, list_cMir, dict_sExpr)
    dict_sTargetPos = get_targetpos_by_gene (dict_sTargets)

    # Open final rbp-hotspot file
    sInDir          = '%s/06_rbp_hotspots/%sbps/%s'   % (sBASE_DIR, nClustDist, nBackWinSize)
    sInFile         = '%s/rbphotspot_%s.dat'          % (sInDir, sFileTag)
    InFile          = open(sInFile, 'rb')
    list_cRBP       = pickle.load(InFile)
    InFile.close()

    for cRBP in list_cRBP:

        cRBP.nOverlap_mir = None # Arbituary flag for None

        try: list_sTargets = dict_sTargetPos[cRBP.sGeneSym]
        except KeyError: continue

        for sPosKey in list_sTargets:
            sChrID       = sPosKey.split(':')[0]
            nStart, nEnd = [int(sPos) for sPos in sPosKey.split(':')[1].split('-')]

            assert cRBP.sChrID == sChrID

            list_nPos_rpb = [cRBP.nStartPos + i for i in range(cRBP.nEndPos - cRBP.nStartPos + 1)]
            list_nPos_mir = [nStart + i for i in range(nEnd - nStart + 1)]
            list_nOverlap = list(set(list_nPos_rpb) & set(list_nPos_mir))

            if list_nOverlap: cRBP.nOverlap_mir = len(list_nOverlap)
            else:
                #With BufferSize
                list_nPos_rpb = [(cRBP.nStartPos - nBufferSize) + i for i in range(cRBP.nEndPos + nBufferSize - cRBP.nStartPos + 1)]
                list_nPos_mir = [nStart + i for i in range(nEnd - nStart + 1)]
                list_nOverlap = list(set(list_nPos_rpb) & set(list_nPos_mir))

                if list_nOverlap: cRBP.nOverlap_mir = 'Buff'
        #loop END: sPosKey
    #loop END: cRBP
    '''
    print(len(list_cRBP))
    dict_sTest = {}
    for cRBP in list_cRBP:

        if not cRBP.nOverlap_mir:
            sTag = cRBP.nOverlap_mir
        elif cRBP.nOverlap_mir == 'Buff':
            sTag = cRBP.nOverlap_mir
        else:
            sTag = 'Mir'

        if sTag not in dict_sTest:
            dict_sTest[sTag] = 0
        dict_sTest[sTag] += 1

    for e in dict_sTest:
        print(e, dict_sTest[e])
    sys.exit()
    '''

    print('Dumped RBPs', len(list_cRBP))

    sOutFile = '%s/rbphotspot_wMir_%s.dat' % (sOutDir, sFileTag)
    OutFile  = open(sOutFile, 'wb')
    pickle.dump(list_cRBP, OutFile)
    OutFile.close()
    OutFile  = open(sOutFile.replace('.dat', '.txt'), 'w')
    for i, cRBP in enumerate(list_cRBP):
        sOut = ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                i, cRBP.nClusterID, cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos,
                cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nPosCnt, cRBP.fQvalue,
                dict_sGENIC[cRBP.sLocale], cRBP.sGeneSym, cRBP.sRBPID.split('_')[0], cRBP.nOverlap, cRBP.nOverlap_mir))
        OutFile.write(sOut)
    #loop END: cRBP
    OutFile.close()
#def END: rbp_enrichment


def output_rbpdata_forbuffer (nClustDist, nBackWinSize, nBufferSize, sFileTag):
    print(green('STEP 9-1: Overlap with Bedtools Part 2', 1))
    sOutDir       = '%s/08_rbps_wInters/%sbps/%s'     % (sBASE_DIR, nClustDist, nBackWinSize)
    os.makedirs(sOutDir, exist_ok=True)

    # Open rbp-hotspot file
    sInDir        = '%s/07_rbps_wMirTargets/%sbps/%s' % (sBASE_DIR, nClustDist, nBackWinSize)
    sInFile       = '%s/rbphotspot_wMir_%s.dat'       % (sInDir, sFileTag)
    print('Load Pickled RBP File', sInFile)
    InFile        = open(sInFile, 'rb')
    list_cRBP     = pickle.load(InFile)
    InFile.close()

    # Determine Buffer overlap
    print('Output for Buffer overlap')
    list_cOverlap    = []
    list_cNoOverlap  = []
    for cRBP in list_cRBP:
        if not cRBP.sRBPID and not cRBP.nOverlap_mir:
            list_cNoOverlap.append(cRBP)
        else:
            list_cOverlap.append(cRBP)
    #loop END: cRBP

    print('list_cOverlap', len(list_cOverlap))
    print('list_cNoOverlap', len(list_cNoOverlap))

    '''
    sOutFile        = '%s/RBP.qvalues.test.txt' % sBASE_DIR
    OutFile         = open(sOutFile, 'w')
    list_fQvalues   = [-math.log10(cRBP.fQvalue + fPSEUDO_COUNT) for cRBP in list_cOverlap]
    for fQvalue in list_fQvalues:
        sOut = '%s\t%s\n' % ('0', fQvalue)
        OutFile.write(sOut)
    list_fQvalues = [-math.log10(cRBP.fQvalue + fPSEUDO_COUNT) for cRBP in list_cNoOverlap]
    for fQvalue in list_fQvalues:
        sOut = '%s\t%s\n' % ('1', fQvalue)
        OutFile.write(sOut)
    OutFile.close()
    sys.exit()
    '''

    list_cNoOverlap = sorted(list_cNoOverlap, key=lambda c:(c.sChrID.replace('chr',''), c.nStartPos))
    sOutFile        = '%s/rbphotspot_forBuffer_%s.txt' % (sOutDir, sFileTag)
    OutFile         = open(sOutFile, 'w')
    for cRBP in list_cNoOverlap:

        sPosList  = ','.join([str(nPos) for nPos in cRBP.list_nPos])
        sOut      = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                   (cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos, cRBP.nClusterID,
                    cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nPosCnt,
                    cRBP.fQvalue, cRBP.fBackProb, cRBP.sLocale, cRBP.sGeneSym,
                    cRBP.sGeneID, sPosList)

        OutFile.write(sOut)
    #loop END: cRBP

    OutFile.close()

    sOutFile = '%s/rbphotspot_overlap_%s.dat' % (sOutDir, sFileTag)
    OutFile  = open(sOutFile, 'wb')
    pickle.dump(list_cOverlap, OutFile)
    OutFile.close()
#def END: output_rbpdata_forbuffer


def consolidate_rbplist (list_cRBP):

    list_cOutput     = []
    dict_nClusterID  = {}
    for cRBP in list_cRBP:
        sKey = cRBP.nClusterID
        if sKey not in dict_nClusterID:
            dict_nClusterID[sKey] = []
        dict_nClusterID[sKey].append(cRBP)
    #loop END: cRBP

    for nClusterID in dict_nClusterID:
        if len(dict_nClusterID[nClusterID]) > 1:

            list_sRBPID       = list(set([cRBP.sRBPID.split('_')[0] for cRBP in dict_nClusterID[nClusterID] if cRBP.sRBPID != '.']))
            list_sRBPPos      = list(set(['%s-%s' %(cRBP.nRBPStart, cRBP.nRBPEnd) for cRBP in dict_nClusterID[nClusterID] if cRBP.sRBPID != '.']))
            list_nOverlap     = list(set(['%s' % cRBP.nOverlap for cRBP in dict_nClusterID[nClusterID] if cRBP.sRBPID != '.']))
            cRBP              = dict_nClusterID[nClusterID][0]
            cRBP.sRBPID       = ','.join(list_sRBPID)
            cRBP.nOverlap     = ','.join(list_nOverlap)
            cRBP.list_sRBPPos = ','.join(list_sRBPPos)
            list_cOutput.append(cRBP)
        else:
            cRBP             = dict_nClusterID[nClusterID][0]
            list_cOutput.append(cRBP)
    #loop END: nClusterID

    list_cOutput = sorted(list_cOutput, key=lambda  c:c.fQvalue)
    '''
    for cRBP in list_cOutput:
        sOut = ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        cRBP.nClusterID, cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos,
        cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nCovered, cRBP.fQvalue,
        dict_sGENIC[cRBP.sLocale], cRBP.sGeneSym, cRBP.sRBPID, cRBP.nOverlap, cRBP.nOverlap_mir))
        print(sOut[:-1])
    '''

    return list_cOutput
#def END: consolidate_rbplist


def qsub_determine_rbp_coverage (list_cClipData, nClustDist, nBackWinSize, nBufferSize, sFileTag, sTempScript, sQueue, bTestRun):
    print(green('STEP 9-3: Mark RBPs that overlap with capture region', 1))

    sOutDir         = '%s/09_rbps_capture_overlap/%sbps/%s' % (sBASE_DIR, nClustDist, nBackWinSize)
    sBufferFile     = '%s/rbphotspot_buffer_%s.dat'         % (sOutDir, sFileTag)
    os.makedirs(sOutDir, exist_ok=True)

    '''Load Buffered RBPs'''
    list_cBuffered  = get_rbp_buffered (list_cClipData, nClustDist, nBackWinSize, nBufferSize, sBufferFile)

    # Open previous rbp-hotspot file with overlapped and mir data
    sInDir          = '%s/08_rbps_wInters/%sbps/%s'  % (sBASE_DIR, nClustDist, nBackWinSize)
    sInFile         = '%s/rbphotspot_overlap_%s.dat' % (sInDir, sFileTag)
    InFile          = open(sInFile, 'rb')
    list_cRBP       = pickle.load(InFile) + list_cBuffered
    InFile.close()
    print('Overlap Number of RBP-Cluster Interactions', len(list_cRBP))

    bCapture        = 1
    list_sCapture   = parse_capture_bedfile(sCAPTUREFILE)
    dict_cRefSeq    = cRef_parse_refflat_line(sREFSEQFILE, list_sCapture, bCapture)
    #dict_sGeneSym   = {dict_cRefSeq[sNMID].sGeneSym:dict_cRefSeq[sNMID] for sNMID in dict_cRefSeq}


    for cRBP in list_cRBP:
        try: cRBP.bBuff
        except AttributeError:cRBP.bBuff = 0
        try: cRBP.nOverlap_mir
        except AttributeError:cRBP.nOverlap_mir = None


    '''

    dict_sTest = {}
    dict_nIntergenic  = {}
    for sNMID in dict_cRefSeq:
        cRef = dict_cRefSeq[sNMID]

        if cRef.sChrID not in dict_nIntergenic:
            dict_nIntergenic[cRef.sChrID] = []
        dict_nIntergenic[cRef.sChrID] += cRef.list_nIntergenic
    #loop END: sNMID

    nCnt = 0
    nCnt1 = 0
    for cRBP in list_cRBP:
        if cRBP.sChrID == 'chrM':  continue

        if dict_sGENIC_alt[cRBP.sLocale] == 'intergenic':
        #if cRBP.sLocale == 'intergenic':
            list_nInterPos = sorted(dict_nIntergenic[cRBP.sChrID], key=lambda e:e[0])
            list_nCoverPos = check_covered_variantpos (list_nInterPos, cRBP.list_nPos)

            if not list_nCoverPos: nCnt += 1
            else:
                cRBP.list_nPos_cap    = list_nCoverPos
                cRBP.list_sRBPPos_cap = check_rbp_covered_pos(list_nInterPos, cRBP.list_sRBPPos)

                print(cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos, cRBP.nClusterID,
                    cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nPosCnt,
                    cRBP.fQvalue, cRBP.fBackProb, cRBP.sLocale, cRBP.sGeneSym,
                    cRBP.sGeneID, cRBP.list_sRBPPos, cRBP.list_sRBPPos_cap)
                nCnt1 += 1

    print(nCnt)
    print(nCnt1)

    nCnt = 0
    dict_sTest = {}
    print(len(dict_cRefSeq))
    for cRBP in list_cRBP:
        try: cRBP.bBuff
        except AttributeError:
            cRBP.bBuff = 0

        if not cRBP.bBuff:
            try: dict_sGeneSym[cRBP.sGeneSym]
            except KeyError:
                nCnt += 1
                if cRBP.sLocale  not in dict_sTest:
                    dict_sTest[cRBP.sLocale] = []
                dict_sTest[cRBP.sLocale].append(cRBP)

    #loop END: cRBP
    print(nCnt)
    for e in dict_sTest:
        print(e, len(dict_sTest[e]), np.median([cRBP.nPosCnt for cRBP in dict_sTest[e]]))
    sys.exit()
    '''

    sPickleFile     = '%s/rbphotspot_forcapture_%s.dat' % (sOutDir, sFileTag)
    OutFile         = open(sPickleFile, 'wb')
    pickle.dump(list_cRBP, OutFile)
    OutFile.close()

    sJobName1       = 'Jinman.Capture.Overlap'
    sCaptureDir     = '%s/capture_%s'            % (sOutDir, sFileTag)
    sLogDir         = '%s/log/%s/%s'             % (sBASE_DIR, sJobName1, sTIME_STAMP)
    os.makedirs(sCaptureDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    nTotalRBPs      = len(list_cRBP)
    nJobs           = 500
    list_sOutFiles  = []
    for nJobIndex in range(nJobs):
        ## MD
        nStartIndex      = int(nTotalRBPs*(nJobIndex+0)/nJobs)
        nEndIndex        = int(nTotalRBPs*(nJobIndex+1)/nJobs) - 1
        sOutFile         = '%s/%04d.dat'               % (sCaptureDir, nJobIndex+1)
        list_sOutFiles.append(sOutFile)
        sScript          = 'date; %s check_capture_v2'    % sTempScript
        sScript         += ' %s %s %s %s; date'        % (sPickleFile, nStartIndex, nEndIndex, sOutFile)
        sScript          = 'echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%04d'\
                            % (sScript, sLogDir, sQueue, sJobName1, nJobIndex+1)
        if bTestRun: print(sScript)
        else:        os.system(sScript)
    #loop ENd: nJobIndex

    sOutFile    = '%s/rbphotspot_capture_%s.dat'       % (sOutDir, sFileTag)
    sJobName2   = 'Jinman.Capture.Overlap.Combine'
    sLogDir     = '%s/log/%s/%s'                       % (sBASE_DIR, sJobName2, sTIME_STAMP)
    os.makedirs(sLogDir, exist_ok=True)
    sScript     = 'date; %s combine_capture'           % sTempScript
    sScript    += ' %s %s %s ; date'                   % (sCaptureDir, nJobs, sOutFile)
    sScript     = 'echo "%s" | qsub -cwd -j y -o %s -q %s -hold_jid %s.* -N %s '\
                % (sScript, sLogDir, sQueue, sJobName1, sJobName2)
    if bTestRun: print(sScript)
    else:        os.system(sScript)
#def END: mark_capture_overlap


def check_capture_v2 (sInFile, nStartIndex, nEndIndex, sOutFile):
    bCapture          = 1
    list_sCapture     = parse_capture_bedfile(sCAPTUREFILE)
    dict_cRefSeq      = cRef_parse_refflat_line(sREFSEQFILE, list_sCapture, bCapture)
    dict_sGeneSym     = {dict_cRefSeq[sNMID].sGeneSym:dict_cRefSeq[sNMID] for sNMID in dict_cRefSeq}
    # Get capture overlap intergenic positions
    dict_sAllPos      = {}
    dict_nIntergenic  = {}
    for sNMID in dict_cRefSeq:
        cRef          = dict_cRefSeq[sNMID]

        list_sAllPos  = cRef.list_nIntergenic + cRef.list_nIntrons
        list_sAllPos  += [[nExons, nExonE] for nExons, nExonE in zip(cRef.list_nExonS, cRef.list_nExonE)]

        if cRef.sChrID not in dict_nIntergenic:
            dict_nIntergenic[cRef.sChrID] = []

        if cRef.sChrID not in dict_sAllPos:
            dict_sAllPos[cRef.sChrID] = []

        dict_nIntergenic[cRef.sChrID] += cRef.list_nIntergenic
        dict_sAllPos[cRef.sChrID]     += list_sAllPos
    #loop END: sNMID

    InFile            = open(sInFile, 'rb')
    list_cRBP         = pickle.load(InFile)[int(nStartIndex):int(nEndIndex)]
    InFile.close()

    list_cOutput      = []
    nCnt = 0
    print('Starting RBPs', len(list_cRBP))

    for cRBP in list_cRBP:

        if cRBP.sChrID == 'chrM':  continue

        bPass = 1
        if dict_sGENIC_alt[cRBP.sLocale] != 'intronic':
            list_nPos      = sorted(dict_sAllPos[cRBP.sChrID], key=lambda e:e[0])
            list_nCoverPos = check_covered_variantpos (list_nPos, cRBP.list_nPos)

            if not list_nCoverPos: bPass = 0
            else:
                cRBP.list_nPos_cap    = list_nCoverPos
                cRBP.list_sRBPPos_cap = check_rbp_covered_pos(list_nPos, cRBP.list_sRBPPos)
        else:
            for sGeneSym in cRBP.sGeneSym.split(','):  #Some genesyms have two genes
                try:             cRef = dict_sGeneSym[sGeneSym]
                except KeyError: bPass = 0

                list_nCoverPos = check_covered_variantpos(cRef.list_nIntrons, cRBP.list_nPos)
                if not list_nCoverPos: bPass = 0
                else:
                    cRBP.list_nPos_cap    = list_nCoverPos
                    cRBP.list_sRBPPos_cap = check_rbp_covered_pos(cRef.list_nIntrons, cRBP.list_sRBPPos)
            #loop END: sGeneSym
        #if END
        if bPass:
            #if len(cRBP.list_sRBPPos.split(',')) != len(cRBP.list_sRBPPos_cap.split(',')):
            #print('************************', cRBP.sGeneSym, dict_sGENIC_alt[cRBP.sLocale], cRBP.bBuff)
            #print('Before', cRBP.list_nPos)
            #print('After', cRBP.list_nPos_cap)
            #print('Before', cRBP.list_sRBPPos,
            #      sum([int(nPosKey.split('-')[1]) - int(nPosKey.split('-')[0]) for nPosKey in cRBP.list_sRBPPos.split(',')]))
            #print('After', cRBP.list_sRBPPos_cap,
            #      sum([int(nPosKey.split('-')[1]) - int(nPosKey.split('-')[0]) for nPosKey in cRBP.list_sRBPPos_cap.split(',')]))
            list_cOutput.append(cRBP)
        else: nCnt += 1
    print('NotFound', nCnt)
    print('Final RBPs', len(list_cOutput))

    OutFile         = open(sOutFile, 'wb')
    pickle.dump(list_cOutput, OutFile)
    OutFile.close()

    OutFile         = open(sOutFile.replace('.dat', '.txt'), 'w')
    for cRBP in list_cOutput:
        sOut      = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                   (cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos, cRBP.nClusterID,
                    cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nPosCnt,
                    cRBP.fQvalue, cRBP.fBackProb, dict_sGENIC_alt[cRBP.sLocale], cRBP.sGeneSym,
                    cRBP.sGeneID, cRBP.list_sRBPPos, cRBP.list_sRBPPos_cap)
        OutFile.write(sOut)
    OutFile.close()
#def END: check_capture_v2


def check_capture (sInFile, nStartIndex, nEndIndex, sOutFile):
    bCapture          = 1
    list_sCapture     = parse_capture_bedfile(sCAPTUREFILE)
    dict_cRefSeq      = cRef_parse_refflat_line(sREFSEQFILE, list_sCapture, bCapture)
    dict_sGeneSym     = {dict_cRefSeq[sNMID].sGeneSym:dict_cRefSeq[sNMID] for sNMID in dict_cRefSeq}
    # Get capture overlap intergenic positions
    dict_sAllPos      = {}
    dict_nIntergenic  = {}
    for sNMID in dict_cRefSeq:
        cRef          = dict_cRefSeq[sNMID]

        list_sAllPos  = cRef.list_nIntergenic + cRef.list_nIntrons
        list_sAllPos  += [[nExons, nExonE] for nExons, nExonE in zip(cRef.list_nExonS, cRef.list_nExonE)]

        if cRef.sChrID not in dict_nIntergenic:
            dict_nIntergenic[cRef.sChrID] = []

        if cRef.sChrID not in dict_sAllPos:
            dict_sAllPos[cRef.sChrID] = []

        dict_nIntergenic[cRef.sChrID] += cRef.list_nIntergenic
        dict_sAllPos[cRef.sChrID]     += list_sAllPos
    #loop END: sNMID

    InFile            = open(sInFile, 'rb')
    list_cRBP         = pickle.load(InFile)[int(nStartIndex):int(nEndIndex)]
    InFile.close()

    list_cOutput      = []
    nCnt = 0
    print('Starting RBPs', len(list_cRBP))
    for cRBP in list_cRBP:

        if cRBP.sChrID == 'chrM':  continue

        bPass = 1
        if dict_sGENIC_alt[cRBP.sLocale] == 'intergenic':
        #if cRBP.sLocale == 'intergenic':
            list_nInterPos = sorted(dict_nIntergenic[cRBP.sChrID], key=lambda e:e[0])
            list_nCoverPos = check_covered_variantpos (list_nInterPos, cRBP.list_nPos)

            if not list_nCoverPos: bPass = 0
            else:
                cRBP.list_nPos_cap    = list_nCoverPos
                cRBP.list_sRBPPos_cap = check_rbp_covered_pos(list_nInterPos, cRBP.list_sRBPPos)
        else:
            for sGeneSym in cRBP.sGeneSym.split(','):  #Some genesyms have two genes
                try:             cRef = dict_sGeneSym[sGeneSym]
                except KeyError: bPass = 0

                if dict_sGENIC_alt[cRBP.sLocale] == 'intronic':
                    list_nCoverPos = check_covered_variantpos(cRef.list_nIntrons, cRBP.list_nPos)
                    if not list_nCoverPos: bPass = 0
                    else:
                        cRBP.list_nPos_cap    = list_nCoverPos
                        cRBP.list_sRBPPos_cap = check_rbp_covered_pos(cRef.list_nIntrons, cRBP.list_sRBPPos)

                else:
                    list_nExon = [[nExons, nExonE] for nExons, nExonE in zip (cRef.list_nExonS, cRef.list_nExonE)]
                    list_nCoverPos = check_covered_variantpos(list_nExon, cRBP.list_nPos)

                    if not list_nCoverPos: bPass = 0
                    else:
                        cRBP.list_nPos_cap    = list_nCoverPos
                        cRBP.list_sRBPPos_cap = check_rbp_covered_pos(list_nExon, cRBP.list_sRBPPos)
                #if END:
            #loop END: sGeneSym
        # if END:
        if bPass:
            #if len(cRBP.list_sRBPPos.split(',')) != len(cRBP.list_sRBPPos_cap.split(',')):
            print('************************',cRBP.sGeneSym, dict_sGENIC_alt[cRBP.sLocale], cRBP.bBuff)
            print('Before', cRBP.list_nPos)
            print('After', cRBP.list_nPos_cap)
            if cRBP.list_sRBPPos:
                print('Before', cRBP.list_sRBPPos, sum([int(nPosKey.split('-')[1]) - int(nPosKey.split('-')[0]) for nPosKey in cRBP.list_sRBPPos.split(',')]))
                print('After', cRBP.list_sRBPPos_cap, sum([int(nPosKey.split('-')[1]) - int(nPosKey.split('-')[0]) for nPosKey in cRBP.list_sRBPPos_cap.split(',')]))
            list_cOutput.append(cRBP)
        else: nCnt += 1
    # loop END: cRBP
    print('NotFound', nCnt)
    print('Final RBPs', len(list_cOutput))

    sys.exit()

    OutFile         = open(sOutFile, 'wb')
    pickle.dump(list_cOutput, OutFile)
    OutFile.close()

    OutFile         = open(sOutFile.replace('.dat', '.txt'), 'w')
    for cRBP in list_cOutput:
        sOut      = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                   (cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos, cRBP.nClusterID,
                    cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nPosCnt,
                    cRBP.fQvalue, cRBP.fBackProb, dict_sGENIC_alt[cRBP.sLocale], cRBP.sGeneSym,
                    cRBP.sGeneID, cRBP.list_sRBPPos, cRBP.list_sRBPPos_cap)
        OutFile.write(sOut)
    OutFile.close()
#def check_capture


def check_covered_variantpos(list_nCover, list_nPos):
    list_sOutput  = []
    for nStart, nEnd in list_nCover:
        for nVarPos in list_nPos:
            if not nStart <= nVarPos <= nEnd: continue
            list_sOutput.append(nVarPos)
        #loop END: nVarPos
    #loop END: nStart, nEnd
    return list(set(list_sOutput))
#def END: get_overlap_length


def check_rbp_covered_pos (list_nCover, list_sRBPPos):

    list_sOutput   = []
    list_sRBPPos   = list_sRBPPos.split(',')
    for nStart, nEnd in list_nCover:

        for sRBPPos in list_sRBPPos:
            if not sRBPPos: continue
            nRBPStart, nRBPEnd = [int(sPos) for sPos in sRBPPos.split('-')]

            if nEnd < nRBPStart: continue
            if nStart > nRBPEnd: continue

            #print('Cover', nStart, nEnd, nEnd - nStart)
            #print('RBP', nRBPStart, nRBPEnd, nRBPEnd- nRBPStart)
            #print('Overlap', max(nStart, nRBPStart), min(nEnd, nRBPEnd), min(nEnd, nRBPEnd) - max(nStart, nRBPStart))

            list_sOutput.append([max(nStart, nRBPStart), min(nEnd, nRBPEnd)])
        #loop END: sRBPPos
    #loop END: sPosKey
    list_sOutput = sorted(list_sOutput, key= lambda e:e[0])

    if not list_sOutput: return ''
    else: return ','.join(['%s-%s' % (nStart, nEnd) for nStart, nEnd in list_sOutput if nStart != nEnd])
#def END: check_rbp_covered_pos

'''
def check_capture_old (sInFile, nStartIndex, nEndIndex, sOutFile):
    bCapture        = 1
    list_sCapture   = parse_capture_bedfile(sCAPTUREFILE)
    dict_cRefSeq    = cRef_parse_refflat_line(sREFSEQFILE, list_sCapture, bCapture)

    InFile          = open(sInFile, 'rb')
    list_cRBP       = pickle.load(InFile)[int(nStartIndex):int(nEndIndex)]
    InFile.close()

    list_cOutput    = []
    print('Starting RBPs', len(list_cRBP))
    for cRBP in list_cRBP:
        if cRBP.sRBPID == '.':       list_cOutput.append(cRBP)
        elif cRBP.sChrID == 'chrM':  list_cOutput.append(cRBP)
        else:
            list_nCoverPos = []
            list_sCapture  = dict_sCapture[cRBP.sChrID]
            list_sRBPPos   = cRBP.list_sRBPPos.split(',')

            for sPosKey in list_sCapture:
                nStart, nEnd  = [int(sPos) for sPos in sPosKey.split(',')]
                for sRBPPos in list_sRBPPos:
                    if not sRBPPos: continue
                    nRBPStart, nRBPEnd = [int(sPos) for sPos in sRBPPos.split('-')]

                    if nEnd < nRBPStart: continue
                    if nStart > nRBPEnd: continue

                    list_nCoverPos.append([max(nStart, nRBPStart), min(nEnd, nRBPEnd)])
                #loop END: sRBPPos
            #loop END: sPosKey
            if not list_nCoverPos: continue

            #print(cRBP.sLocale, len(list_nCoverPos), cRBP.sRBPID)
            list_nVarCover = []
            for nStart, nEnd in list_nCoverPos:

                if nEnd < cRBP.nStartPos: continue
                if nStart > cRBP.nEndPos: continue

                #if cRBP.nStartPos >= nStart and cRBP.nEndPos <= nEnd:
                list_nVarCover.append([nStart,nEnd])
            #loop END: nStart, nEnd
            if list_nVarCover:
                cRBP.list_sRBPPos_cap = ','.join(['%s-%s' % (nStart, nEnd) for nStart, nEnd in list_nVarCover])
                list_cOutput.append(cRBP)
        #if END:
    #loop END: cRBP

    print('Final RBPs', len(list_cOutput))

    OutFile         = open(sOutFile, 'wb')
    pickle.dump(list_cOutput, OutFile)
    OutFile.close()
#def check_capture
'''


def combine_capture (sCaptureDir, nJobs, sOutFile):
    list_cRBP = []
    for nJobIndex in range(int(nJobs)):
        ## MD
        sInFile       = '%s/%04d.dat' % (sCaptureDir, nJobIndex+1)
        InFile        = open(sInFile, 'rb')
        list_cRBP    += pickle.load(InFile)
        InFile.close()
    #loop END: nJobIndex
    OutFile = open(sOutFile, 'wb')
    pickle.dump(list_cRBP, OutFile)
    OutFile.close()
#def END: combine_capture


def determine_rbp_coverage (list_cClipData, nClustDist, nBackWinSize, nBufferSize, sFileTag):
    print(green('STEP 10: Determine RBP coverage ', 1))

    bCapture          = 1
    list_sCapture     = parse_capture_bedfile(sCAPTUREFILE)
    dict_cRefSeq      = cRef_parse_refflat_line(sREFSEQFILE, list_sCapture, bCapture)
    dict_sGenicSize   = calculate_genic_size (dict_cRefSeq, bCapture)

    # dict_sGeneSym     = {dict_cRefSeq[sNMID].sGeneSym:dict_cRefSeq[sNMID] for sNMID in dict_cRefSeq}

    # Open previous rbp-hotspot file with overlapped and mir data
    sInDir          = '%s/09_rbps_capture_overlap/%sbps/%s'  % (sBASE_DIR, nClustDist, nBackWinSize)
    sInFile         = '%s/rbphotspot_capture_%s.dat'   % (sInDir, sFileTag)
    InFile          = open(sInFile, 'rb')
    list_cRBP       = pickle.load(InFile)
    InFile.close()
    print('Number of RBP-Cluster Interactions', len(list_cRBP))

    dict_sLocale    = {}
    for cRBP in list_cRBP:
        sKey = dict_sGENIC_alt[cRBP.sLocale]
        if sKey not in dict_sLocale:
            dict_sLocale[sKey] = []
        if bCapture:
            dict_sLocale[sKey] += cRBP.list_sRBPPos_cap.split(',')
        else:
            dict_sLocale[sKey] += cRBP.list_sRBPPos.split(',')
    #loop END: cRBP

    for sLocale in dict_sLocale:
        try: nGenicLen    = dict_sGenicSize[sLocale]
        except KeyError: continue #Skips 'others'

        list_sRBPLen = list(filter(None, list(set(dict_sLocale[sLocale]))))
        nTotRBPLen   = sum([int(nPosKey.split('-')[1]) - int(nPosKey.split('-')[0]) + 1 for nPosKey in list_sRBPLen])

        #if sLocale == 'intergenic':
            #print(sLocale, nTotRBPLen, nGenicLen)
        #else:
        print(sLocale, nTotRBPLen, nGenicLen, nTotRBPLen/nGenicLen)

    '''print(cRBP.nClusterID, cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos,
          cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nPosCnt, cRBP.fQvalue,
          cRBP.sLocale, cRBP.sGeneSym, cRBP.sRBPID, cRBP.nRBPStart, cRBP.nRBPEnd, cRBP.list_sRBPPos,
          cRBP.nOverlap, cRBP.nOverlap_mir)
    '''
#def END: mark_capture_overlap


def normalize_enrichment (nClustDist, nBackWinSize, nBufferSize, fRecurrency, sFileTag):
    print(green('STEP 11: Determine Enrichment by Genic Region', 1))
    bCapture        = 1
    fPvalue         = 0.05
    bRecurrency     = 1
    list_sCapture   = parse_capture_bedfile(sCAPTUREFILE)
    dict_cRefSeq    = cRef_parse_refflat_line(sREFSEQFILE, list_sCapture, bCapture)
    dict_sGenicSize = calculate_genic_size (dict_cRefSeq, bCapture)

    # Open previous rbp-hotspot file with overlapped and mir data
    sInDir          = '%s/09_rbps_capture_overlap/%sbps/%s'  % (sBASE_DIR, nClustDist, nBackWinSize)
    sInFile         = '%s/rbphotspot_capture_%s.dat'         % (sInDir, sFileTag)
    InFile          = open(sInFile, 'rb')
    list_cRBP       = pickle.load(InFile)
    InFile.close()
    print('Number of RBP-Cluster Interactions', len(list_cRBP))

    print('Classify by Overlap')
    dict_sLocale    = {}

    for cRBP in list_cRBP:
        if fPvalue:
            if cRBP.fQvalue > fPvalue: continue

        if bRecurrency:
            if cRBP.nPatCnt <= int(nSAMPLE_CNT * fRecurrency): continue


        if cRBP.bBuff:
            if not cRBP.sRBPID: bRBP = 'NoRBP'
            else: bRBP = 'Buff'
        else: bRBP = 'RBP'

        if not cRBP.nOverlap_mir:
            bMir = 'NoMir'
        else:
            try:
                bMir = 'Mir' if cRBP.nOverlap_mir > 0 else cRBP.nOverlap_mir
            except TypeError:
                bMir = cRBP.nOverlap_mir

        cRBP.sCategory  = '%s-%s' % (bRBP, bMir)
        sKey            = dict_sGENIC_alt[cRBP.sLocale]
        if sKey not in dict_sLocale:
            dict_sLocale[sKey] = []
        dict_sLocale[sKey].append(cRBP)
    # loop END: cRBP


    dict_sCategories = {'NoRBP-NoMir':'0', 'RBP-NoMir':'1', 'NoRBP-Mir': '2',
                        'RBP-Mir':'3', 'NoRBP-Buff':'4', 'Buff-NoMir':'5','RBP-Buff':'6'}
    dict_sOutputKey  = {'NoRBP-NoMir':'NoRBP-NoMir', 'RBP-NoMir':'RBP-NoMir', 'NoRBP-Mir': 'NoRBP-Mir',
                        'RBP-Mir':'RBP-Mir', 'NoRBP-Buff':'Others', 'Buff-NoMir':'Others','RBP-Buff':'Others'}
    list_sClasses    = list(dict_sCategories.keys())
    list_sOutputKey  = ['RBP-NoMir', 'NoRBP-NoMir', 'NoRBP-Mir', 'RBP-Mir', 'Others']
    list_sRegion     = ['intergenic', 'intronic', 'exonic', '5UTR', '3UTR']
    dict_nRBPLen     = get_rbp_cover(list_cRBP, bCapture)


    sHeader          = '\t'.join([sClass for sClass in list_sOutputKey])
    print('region\t%s' % sHeader)
    for sLocale in list_sRegion:

        nGenicLen   = dict_sGenicSize[sLocale]
        nRBPLen     = dict_nRBPLen[sLocale]

        dict_sRBPClass = {}
        for cRBP in dict_sLocale[sLocale]:
            sKey = cRBP.sCategory
            if sKey not in dict_sRBPClass:
                dict_sRBPClass[sKey] = []
            dict_sRBPClass[sKey].append(cRBP)
        #loop END: cRBP

        dict_sOutput = {sClass: 0.0 for sClass in list_sOutputKey}
        for sClass in list_sClasses:

            if dict_sCategories[sClass] in ['0', '2', '4', '5', '6']:
                nDenom = nGenicLen
            else:
                nDenom = nRBPLen

            try:
                list_cRBP    = dict_sRBPClass[sClass]
                ########################
                nCategoryCnt = sum([cRBP.nPosCnt for cRBP in list_cRBP])
                ########################
                if nDenom:
                    nNormCnt     = (nCategoryCnt / nDenom) * 1000000
                else: nNormCnt   = 0

                dict_sOutput[dict_sOutputKey[sClass]] = nDenom

            except KeyError:
                dict_sOutput[dict_sOutputKey[sClass]] = 0.0

        #loop END: sClass

        sOut    = '\t'.join(['%0.4f' % dict_sOutput[sClass] for sClass in list_sOutputKey])
        print('%s\t%s' % (sLocale if sLocale != 'exonic' else 'cds', sOut))


    #loop END: sLocale
# def END: determine_intermediates


def get_rbp_cover (list_cRBP, bCapture):

    dict_sLocale = {}
    for cRBP in list_cRBP:

        sKey = dict_sGENIC_alt[cRBP.sLocale]
        if sKey not in dict_sLocale:
            dict_sLocale[sKey] = []

        if bCapture:
            dict_sLocale[sKey] += cRBP.list_sRBPPos_cap.split(',')
        else:
            dict_sLocale[sKey] += cRBP.list_sRBPPos.split(',')
    #loop END: cRBP
    dict_sOutput = {}
    for sLocale in dict_sLocale:

        list_sRBPLen = list(filter(None, list(set(dict_sLocale[sLocale]))))
        nTotRBPLen   = sum([int(nPosKey.split('-')[1]) - int(nPosKey.split('-')[0]) + 1for nPosKey in list_sRBPLen])

        if sLocale not in dict_sOutput:
            dict_sOutput[sLocale] = 0.0
        dict_sOutput[sLocale] = nTotRBPLen
    #loop END: sLocale

    for e in dict_sOutput:
        print(e, dict_sOutput[e])

    #VS-Check
    if not dict_sOutput:
        sys.exit('Invalid Dict : get_rbp_cover : dict_sOutput size= %d' % len(dict_sOutput))
    return dict_sOutput
#def END: get_rbp_cover


def get_rbp_buffered (list_cClipData, nClustDist, nBackWinSize, nBufferSize, sBufferFile):
    if not os.path.isfile(sBufferFile):

        print('Loading Buffered RBPs', sBufferFile)

        sInDir        = '%s/08_rbps_wInters/%sbps/%s/out%s'% (sBASE_DIR, nClustDist, nBackWinSize, nBufferSize)
        list_cRBP_buf = []
        for cClipData in list_cClipData:
            sInFile         = '%s/%s.intersect.dat' % (sInDir, cClipData.sAccID)
            #list_sIntersect = parse_intersect_file(sInFile)
            InFile          = open(sInFile, 'rb')
            list_sIntersect = pickle.load(InFile)
            InFile.close()
            list_cRBP_buf   += list_sIntersect
        #loop END: cClipData

        for cRBP in list_cRBP_buf:
            cRBP.bBuff = 1  # Arbituary Buffered Flag
        #loop END: cRBP

        list_cBuffered = consolidate_rbplist(list_cRBP_buf)

        print('Compiled buffer list', len(list_cRBP_buf))
        print('Consolidated buffer list', len(list_cBuffered))

        OutFile  = open(sBufferFile, 'wb')
        pickle.dump(list_cBuffered, OutFile)
        OutFile.close()
    else:
        print('Loading Pickle Buffered RBPs', sBufferFile)
        InFile        = open(sBufferFile, 'rb')
        list_cBuffered = pickle.load(InFile)
        InFile.close()
    #if END:
    return list_cBuffered
#def END: get_rbp_buffered


def rbp_enrichment (dict_c3UTRs, nClustDist, nBackWinSize, fRecurrency):
    # Load TCGA FPKM File List
    sFPKMDir        = '%s/%s/FPKMfiles/UQ_norm_recode'    % (sBASE_DIR, sSTUDY)
    sFPKMFile       = '%s/%s/%s_rnaseq-UQ.txt'            % (sBASE_DIR, sSTUDY, sSTUDY)
    dict_sFPKMFiles = load_fpkm_list (sFPKMFile)
    dict_sExpr      = get_expr_data (sFPKMDir, dict_sFPKMFiles)

    # Load Top miRNAs in HepG2
    list_cMir       = load_top_mirnas (dict_sMIRFILE[sCELLS])

    # Open final rbp-hotspot file
    sInDir        = '%s/06_final_rbp_hotspots/%sbps/%s' % (sBASE_DIR, nClustDist, nBackWinSize)
    sInFile       = '%s/final_rbphotspot_wNonOverlap.dat'           % sInDir
    InFile        = open(sInFile, 'rb')
    list_cRBP     = pickle.load(InFile)
    InFile.close()
    '''
    list_cMir     = cMirBaseData_read_GFF3_file (dict_sMIRBASE['hg19'])
    dict_cMir     = {}
    for cMir in list_cMir:
        if cMir.sMirName not in dict_cMir:
            dict_cMir[cMir.sMirName] = ''
        dict_cMir[cMir.sMirName] = cMir.sMirSeq

    for sMirName in dict_cMir:
        if 'mir-128' in sMirName:
            print(sMirName, dict_cMir[sMirName].replace('T', 'U'))

    sys.exit()
    '''



    print(len(list_cRBP))
    list_cOverlap    = []
    list_cNonOverlap = []

    for cRBP in list_cRBP:
        if cRBP.sLocale == 'exonic': continue # Filter exonic
        if fRecurrency > 0:                   # Recurrency Filter ON
            if cRBP.nPatCnt <= int(nSAMPLE_CNT * fRecurrency): continue

        if cRBP.nOverlap == 0: list_cNonOverlap.append(cRBP)
        else:                  list_cOverlap.append(cRBP)
    #loop END: cRBP

    print(len(list_cNonOverlap))
    print(len(list_cOverlap))

    list_cOverlap         = sorted(list_cOverlap,    key=lambda c:c.fQvalue)[:50]
    list_cNonOverlap      = sorted(list_cNonOverlap, key=lambda c:c.fQvalue)[:50]

    list_cNonOverlap_3UTR = [cRBP for cRBP in list_cNonOverlap if cRBP.sLocale == '3UTR']

    print(len(list_cOverlap))
    sys.exit()

    for cRBP in list_cNonOverlap:
        sOut = ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                cRBP.nClusterID, cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos,
                cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nCovered, cRBP.fQvalue,
                dict_sGENIC[cRBP.sLocale], cRBP.sGeneSym, cRBP.sRBPID.split('_')[0], cRBP.nOverlap))

        print(sOut[:-1])

        sScript = 'bedtools intersect -a %s -b %s -wao' % (sHotspotFile, sBedFile)

        #print(-math.log(cRBP.fQvalue +fPSEUDO_COUNT, 10))
#def END: rbp_enrichment










'''
def basic_stats_clusters_rbps_old (list_cClipData, nClustDist, nBackWinSize, fRecurrency):
    print('STEP 7: Basic Statistics of Clusters-RBPs P<0.05')

    ## Load COSMIC data ##
    dict_cCosGene = load_cosmic()
    ## Load CPDB data ##
    dict_sCPDBGene = load_CPDB()

    list_sLocale  = ['exonic', '3UTR', '5UTR', 'intergenic', 'intronic', 'splicing', 'upstream', 'downstream']
    sInDir        = '%s/05_rbp_intersect/%sbps/%s/%s' % (sBASE_DIR, nClustDist, nBackWinSize, 0)
    list_sAllRBPs = []
    for cClipData in list_cClipData:
        sInFile         = '%s/%s.intersect.txt' % (sInDir, cClipData.sAccID)
        list_sIntersect = parse_intersect_file(sInFile)
        list_sAllRBPs  += list_sIntersect
    #loop END: cClipData

    nRBP_HotspotCnt = 0
    dict_sRBP       = {}
    nCnt            = 0
    for cRBP in list_sAllRBPs:

        if fRecurrency > 0: #Recurrency Filter ON
            if cRBP.nPatCnt <= int(nSAMPLE_CNT * fRecurrency): continue

        if cRBP.sGeneSym not in dict_sCPDBGene: #CPDB ovelap
            nCnt += 1
            continue

        if cRBP.sLocale == 'exonic': continue #Filter exonic

        nRBP_HotspotCnt += 1

        sKey = cRBP.sRBPID.split('_')[0]

        if sKey not in dict_sRBP:
            dict_sRBP[sKey] = []
        dict_sRBP[sKey].append(cRBP)

    #loop END: cRBP
    print('Filtered', nCnt)
    print('Number of RBP-Cluster Interactions', nRBP_HotspotCnt)

    for sRBPID in dict_sRBP:
        if len(dict_sRBP[sRBPID]) == 1: continue

        print(sRBPID, len(dict_sRBP[sRBPID]))

        for cRBP in dict_sRBP[sRBPID]:
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (
            cRBP.nClusterID, cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos,
            cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nCovered, cRBP.fQvalue,
            cRBP.sLocale, cRBP.sGeneSym, cRBP.sRBPID, cRBP.nOverlap))

    sys.exit()


    dict_sLocale    = {}
    list_sOutput    = []
    for sRBPID in dict_sRBP:
        dict_nClusterID = {}
        if len(dict_sRBP[sRBPID]) == 1: continue

        for cRBP in dict_sRBP[sRBPID]:
            if cRBP.nClusterID not in dict_nClusterID:
                dict_nClusterID[cRBP.nClusterID] = []
            dict_nClusterID[cRBP.nClusterID].append(cRBP)
        #loop END: cRBP

        for nClusterID in dict_nClusterID:
            if len(dict_nClusterID[nClusterID]) == 1: continue
            #print(nClusterID, len(dict_nClusterID[nClusterID]))

            for cRBP in dict_nClusterID[nClusterID][:1]:

                list_sOutput.append(cRBP)
                sKey = dict_sGENIC[cRBP.sLocale]

                if sKey not in dict_sLocale:
                    dict_sLocale[sKey] = []
                dict_sLocale[sKey].append(cRBP)

            #loop END: cRBP
        #loop END: nClusterID
    #loop END: sRBPID
    print('list_sOutput', len(list_sOutput))
    for cRBP in list_sOutput:
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (cRBP.nClusterID, cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos,
                                                                        cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.nCovered, cRBP.fQvalue,
                                                                        cRBP.sLocale, cRBP.sGeneSym, cRBP.sRBPID, cRBP.nOverlap))


    for sLocale in dict_sLocale:
        print(sLocale, len(dict_sLocale[sLocale]))
        #for cRBP in dict_sLocale[sLocale]:
            #print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (cRBP.sChrID, cRBP.nStartPos, cRBP.nEndPos,
            #                                                        cRBP.nSize, cRBP.nVarCnt, cRBP.nPatCnt, cRBP.fQvalue,
            #                                                        cRBP.sLocale, cRBP.sGeneSym, cRBP.sRBPID, cRBP.nOverlap))

    #loop END: sLocale

    sys.exit()
    #Locale by RBP
    for sRBPName in dict_sRBP:
        dict_sLocale = {sLocale:0 for sLocale in list_sLocale}

        for cRBP in dict_sRBP[sRBPName]:
            sKey = dict_sGENIC[cRBP.sLocale]

            if sKey not in dict_sLocale:
                dict_sLocale[sKey] = 0
            dict_sLocale[sKey] += 1
        #loop END: cRBP
        print('%s\t%s\t%s' % (sRBPName, len(dict_sRBP[sRBPName]),  '\t'.join(str(dict_sLocale[sLocale]) for sLocale in list_sLocale)))
    #loop END: sRBPName
#def END: basic_stats_clusters_rbps


def OLD_qsub_get_rbp_variants (list_cClipData, sTempScript,  sQueue, bTestRun):

    sJobName      = 'Jinman.RBP.variants'
    sOutDir       = '%s/02_rbp_vardata/%s'  % (sBASE_DIR, sSTUDY)
    sLogDir       = '%s/log/%s/%s'          % (sBASE_DIR, sJobName, sTIME_STAMP)

    os.makedirs(sOutDir, exist_ok=True)

    for cClipData in list_cClipData:

        sScript = '%s get_rbp_variants %s %s' % (sTempScript, sOutDir, cClipData.sAccID)
        sScript = 'echo "%s" | qsub -cwd -j y -o %s -q %s -N %s.%s' \
                    % (sScript, sLogDir, sQueue, sJobName, cClipData.sAccID)

        if bTestRun: print(sScript)
        else:
            os.makedirs(sLogDir, exist_ok=True)
            os.system(sScript)
    # loop END: cClipData
#def END: qsub_get_rbp_variants


def OLD_get_rbp_variants (sOutDir, sAccID):

    print_start_msg('Annotating BedFile', sBASE_DIR)

    sClipBedDir    = '%s/sorted_bed'    % sECLIP_DIR
    sBedFile       = '%s/%s.bed'        % (sClipBedDir, sAccID)
    if not os.path.isfile(sBedFile): sys.exit('File Not Found %s' % sBedFile)

    sVCFDir        = '%s/%s/VCFs_hg19'  % (sBASE_DIR, sSTUDY)

    # Load Patient IDs
    list_sPatIDs   = get_patientIDs(sSTUDY, sSEQDATA)

    dict_cVCF      = {}
    for sPatID in list_sPatIDs:

        sVCFFile   = '%s/%s.vcf.gz'   % (sVCFDir, sPatID)
        sScript    = 'tabix %s -R %s' % (sVCFFile, sBedFile)

        sStdOut    = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout
        for i, sReadLine in enumerate(sStdOut):

            sVarInfo            = str(sReadLine, 'UTF-8').strip('\n')
            list_sColumn        = sVarInfo.split('\t')

            cVCF                = cVCFData()
            cVCF.sPatID         = sPatID
            cVCF.sChrID         = list_sColumn[0]
            cVCF.nPos           = int(list_sColumn[1])
            cVCF.sDBSNP_ID      = list_sColumn[2]
            cVCF.sRefNuc        = list_sColumn[3]
            cVCF.sAltNuc        = list_sColumn[4]
            cVCF.fQual          = list_sColumn[5]
            cVCF.sFilter        = list_sColumn[6]
            cVCF.sInfo          = list_sColumn[7]
            cVCF.sFormat        = list_sColumn[8]
            cVCF.list_sMisc     = list_sColumn[9:]

            if cVCF.sChrID not in dict_cVCF:
                dict_cVCF[cVCF.sChrID] = []

            dict_cVCF[cVCF.sChrID].append(cVCF)
        #loop END: i, sReadLine
    #loop END: sPatID

    #VS Check
    if not dict_cVCF:
        sys.exit('ERROR: Invalid dict_sVariants Size= %d' % (len(dict_cVCF)))

    #Pickle dump
    OutFile = open('%s/%s.var' % (sOutDir, sAccID), 'wb')
    pickle.dump(dict_cVCF, OutFile)
    OutFile.close()

    print_done_msg('Annotating BedFile', sBASE_DIR)
#def END: get_rbp_variants


def filter_rbp_sites (list_cClipData, dict_cRefSeq):
    sVarDir     = '%s/02_rbp_vardata/%s'    % (sBASE_DIR, sSTUDY)
    sClipBedDir = '%s/sorted_bed'           % sECLIP_DIR

    for cClipData in list_cClipData[:1]:

        cClipData.sAccID = 'ENCFF025QVR'

        sBedFile         = '%s/%s.bed'        % (sClipBedDir, cClipData.sAccID)
        list_cPeakData   = cPeakData_parse_bedfile(sBedFile, dict_cRefSeq)
        print(len(list_cPeakData))

        sVarFile         = '%s/%s.var'        % (sVarDir, cClipData.sAccID)

        InFile           = open(sVarFile, 'rb')
        dict_cVCF        = pickle.load(InFile)
        InFile.close()

        list_sOutput     = []
        for cPeak in list_cPeakData:

            if cPeak.sChrID not in dict_cVCF: continue

            list_cVCF = dict_cVCF[cPeak.sChrID]

            for cVCF in list_cVCF:

                if cPeak.nStartPos <= cVCF.nPos <= cPeak.nEndPos:
                    cPeak.list_cVCF.append(cVCF)
            #loop END: cVCF

            if not cPeak.list_cVCF: continue
            list_sOutput.append(cPeak)
        #loop END: cPeak
        print(len(list_sOutput))

        for cPeak in sorted(list_sOutput, key=lambda c:len(c.list_cVCF), reverse=True):
            print(cPeak.nStartPos, cPeak.nEndPos, cPeak.nEndPos-cPeak.nStartPos, cPeak.sSample, len(cPeak.list_cVCF))

            for cVCF in cPeak.list_cVCF:
                print(cVCF.sPatID, cVCF.nPos)
            sys.exit()

    #loop END: cCLipData
#def END: filter_rbp_sites

def get_mirseq ():
    list_cMir     = cMirBaseData_read_GFF3_file (dict_sMIRBASE['hg19'])
    dict_cMir     = {}
    for cMir in list_cMir:
        if cMir.sMirName not in dict_cMir:
            dict_cMir[cMir.sMirName] = ''
        dict_cMir[cMir.sMirName] = cMir.sMirSeq

    for sMirName in dict_cMir:
        if 'mir-128' in sMirName:
            print(sMirName, dict_cMir[sMirName].replace('T', 'U'))

    sys.exit()
#def END: get_mirseq

'''





if __name__ == '__main__':
    if len(sys.argv) == 1: main()
    else:
        function_name       = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys(): locals()[function_name](*function_parameters)
        else: sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    #if END: len(sys.argv)
#if END: __name__

