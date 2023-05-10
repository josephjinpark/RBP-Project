#!/extdata6/Jinman/opt/python3/bin/python3

import os, sys, time
import numpy as np

sBASE_DIR   = '/extdata6/Jinman/07_rbp'
sTIME_STAMP = '%s' % (time.ctime().replace(' ', '-').replace(':', '_'))

class cDenovoData: pass

class cRefSeq: pass


def load_cRefSeq (sInFile):
    dict_sOutput = {}
    InFile = open(sInFile, 'r')
    for sReadLine in InFile:
        list_sColumn = sReadLine.strip('\n').split('\t')
        cRef = cRefSeq()
        cRef.sGeneSym = list_sColumn[0].upper()
        cRef.sNMID = list_sColumn[1]

        if cRef.sNMID not in dict_sOutput:
            dict_sOutput[cRef.sGeneSym] = ''
        dict_sOutput[cRef.sGeneSym] = cRef.sNMID
    #loop END: sReadLine

    return dict_sOutput
#def END: load_cRefSeq

## region Denovo-DB

def denovo_db ():
    sWorkDir       = '%s/ref'           % sBASE_DIR
    sInFile        = '%s/denovodb.tsv'  % sWorkDir
    sRefSeqFile    = '/extdata6/Jinman/reference_genome/091613_4WT.txt'
    dict_sRefSeq   = load_cRefSeq(sRefSeqFile)

    dict_sLocalKey = {'3-prime-UTR': '3UTR', 'non-coding-exon': '3UTR',
                      '5-prime-UTR': '5UTR', 'non-coding-exon-near-splice': '5UTR',
                      'missense': 'CDS-nonsyn', 'synonymous': 'CDS-syn', 'coding':'CDS-syn',
                      'intergenic': 'intergenic', 'none': 'intergenic', 'upstream-gene': 'intergenic',
                      'downstream-gene': 'intergenic', 'frameshift': 'intergenic', '': 'intergenic',
                      'intron': 'intronic', 'intron-near-splice': 'intronic',
                      'splice-acceptor': 'intronic', 'splice-donor': 'intronic',
                      'codingComplex': 'intronic', 'codingComplex-near-splice': 'intronic',
                      'coding-near-splice': 'intronic', 'coding-unknown': 'intronic',
                      'missense-near-splice': 'intronic', 'start-lost': 'intronic',
                      'stop-gained': 'intronic', 'stop-gained-near-splice': 'intronic',
                      'stop-lost': 'intronic', 'synonymous-near-splice': 'intronic',
                      'frameshift-near-splice':'intronic'}
    list_sLocale   = ['3UTR', '5UTR', 'CDS-nonsyn', 'CDS-syn', 'intergenic', 'intronic']
    #list_sLocale   = [ 'CDS', '5UTR', '3UTR','intronic', 'intergenic']


    dict_nRefSizes = {'intergenic': 1947984673, 'intronic': 992967918, '5UTR': 3835217,
                      '3UTR': 24213897, 'CDS-nonsyn': 30998295, 'CDS-syn': 30998295}

    dict_sSeqCap2  = {'intergenic': 410433, 'intronic': 18727225, '5UTR': 1278892,
                      '3UTR': 4997154, 'CDS-nonsyn': 26865812, 'CDS-syn': 26865812}

    dict_sAgilent5 = {'intergenic': 146586, 'intronic': 19330436, '5UTR': 756366,
                      '3UTR': 1007900, 'CDS-nonsyn': 27281728, 'CDS-syn': 27281728}

    dict_sAgilent4 = {'intergenic': 333981, 'intronic': 23063098, '5UTR': 476897,
                      '3UTR': 1012396, 'CDS-nonsyn': 23295714, 'CDS-syn': 23295714}

    dict_sAgilent50 = {'intergenic': 154505, 'intronic': 13731775, '5UTR': 710990,
                      '3UTR': 941739, 'CDS-nonsyn': 29366744, 'CDS-syn': 29366744}

    dict_sCapKit   = {'autism':dict_sSeqCap2,
                      'developmental_disorder':dict_sSeqCap2,
                      'congenital_heart_disease':dict_sSeqCap2,
                      'intellectualDisability':dict_sAgilent4,
                      'schizophrenia':dict_sAgilent4,
                      'epilepsy':dict_sSeqCap2}


    list_cDB     = parse_db (sInFile)

    #convert_to_vcf (sWorkDir, list_cDB)
    #annovar_vcf (sWorkDir)
    #sInFile       = '%s/denovodb.anno.vcf.hg19_multianno.vcf' % sWorkDir
    #list_cDB_anno = parse_db_annovar (sInFile)

    list_cDB       = [cDB for cDB in list_cDB if cDB.sPrimaryPhenotype != 'control']
    list_cDB       = [cDB for cDB in list_cDB if cDB.sSequenceType == 'genome']
    print('Total Variants', len(list_cDB))

    #All Variants
    dict_sSample  = {}
    dict_sGeneSym = {}
    dict_sLocale  = {}
    dict_sDisease = {}
    for cDB in list_cDB:

        if cDB.sSampleID not in dict_sSample:
            dict_sSample[cDB.sSampleID]  = []
        dict_sSample[cDB.sSampleID].append(dict_sLocalKey[cDB.sFunctionClass])
        dict_sSample[cDB.sSampleID] = sorted(list(set(dict_sSample[cDB.sSampleID])))

        ### Gene Distribution ###
        sKey = cDB.sGeneSym
        if sKey not in dict_sGeneSym:
           dict_sGeneSym[sKey] = []
        dict_sGeneSym[sKey].append(dict_sLocalKey[cDB.sFunctionClass])
        dict_sGeneSym[sKey] = sorted(list(set(dict_sGeneSym[sKey])))


        ### Variant Distribution ###
        sKey = dict_sLocalKey[cDB.sFunctionClass]
        if sKey not in dict_sLocale:
            dict_sLocale[sKey] = 0
        dict_sLocale[sKey] += 1


        ### Sequence Data Distribution ###
        sKey  = cDB.sPrimaryPhenotype
        sKey2 = cDB.sSequenceType
        sKey3 = cDB.sSampleID
        if sKey not in dict_sDisease:
            dict_sDisease[sKey] = {}

        if sKey2 not in dict_sDisease[sKey]:
            dict_sDisease[sKey][sKey2] = {}

        if sKey3 not in dict_sDisease[sKey][sKey2]:
            dict_sDisease[sKey][sKey2][sKey3] = 0
        dict_sDisease[sKey][sKey2][sKey3] += 1
    #loop END: cDB
    print('Sample Count', len(dict_sSample))
    nCnt = 0
    for sSampleID in dict_sSample:

        if 'CDS-nonsyn' in dict_sSample[sSampleID]: continue
        nCnt += 1
        print(sSampleID, dict_sSample[sSampleID])
    print('No non-syn', nCnt)



    for sDisease in dict_sDisease:
        print(sDisease, len(dict_sDisease[sDisease]), ' '.join(['%s-%s' % (sSeqData, len(dict_sDisease[sDisease][sSeqData])) for sSeqData in dict_sDisease[sDisease]]))
    sys.exit()
    for sLocale in list_sLocale:
        print(sLocale, dict_sLocale[sLocale])

    print(len(dict_sGeneSym))
    dict_sVarDist = {}
    for sGeneSym in dict_sGeneSym:

        sKey = '-'.join(dict_sGeneSym[sGeneSym])

        if sKey not in dict_sVarDist:
            dict_sVarDist[sKey] = []
        dict_sVarDist[sKey].append(sGeneSym)
    #loop END: sGeneSym

    for sKey in dict_sVarDist:
        print(sKey, len(dict_sVarDist[sKey]))


    sys.exit()
    #Disease by Disease
    dict_sDisease = {}
    for cDB in list_cDB:
        sKey = cDB.sPrimaryPhenotype
        if sKey not in dict_sDisease:
            dict_sDisease[sKey] = []
        dict_sDisease[sKey].append(cDB)
    #loop END: cDB

    list_sDisease = ['autism', 'developmental_disorder', 'congenital_heart_disease',
                     'intellectualDisability', 'schizophrenia', 'epilepsy']

    for sDisease in list_sDisease:
        print('%s' % sDisease)
        list_cDB        = dict_sDisease[sDisease]
        dict_sCapLen    = dict_sCapKit[sDisease]

        dict_sLocale    = {}
        for cDB in list_cDB:
            sKey = dict_sLocalKey[cDB.sFunctionClass]
            if sKey not in dict_sLocale:
                dict_sLocale[sKey] = 0
            dict_sLocale[sKey] += 1
        #loop END: cDB
        for sLocale in list_sLocale:
            print(sLocale, dict_sLocale[sLocale], dict_sCapLen[sLocale], (dict_sLocale[sLocale]/dict_sCapLen[sLocale])*1000000 )
    #loop END: sDisease

#def END: main


def parse_db (sInFile):


    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:
        ''' File Format
        0	SampleID	        1-05366
        1	StudyName	        Homsy2015
        2	PubmedID	        26785492
        3	NumProbands	        1213
        4	NumControls	        0
        5	SequenceType	    exome
        6	PrimaryPhenotype	congenital_heart_disease
        7	Validation	        unknown
        8	Chr	                1
        9	Position	        878640
        10	Variant	            A>AG
        11	rsID	            0
        12	DbsnpBuild	        0
        13	AncestralAllele	    N
        14	1000GenomeCount	    0
        15	ExacFreq	        0
        16	EspAaFreq	        0
        17	EspEaFreq	        0
        18	Transcript	        NM_152486.2
        19	codingDnaSize	    2046
        20	Gene	            SAMD11
        21	FunctionClass	    frameshift
        22	cDnaVariant	        c.1577dup
        23	ProteinVariant	    p.L527Tfs*12
        24	Exon/Intron	        exon12
        25	PolyPhen(HDiv)	    -1
        26	PolyPhen(HVar)	    -1
        27	SiftScore	        -1
        28	CaddScore	        -1
        29	LofScore	        1
        30	LrtScore	        -1
        '''
        if sReadLine.startswith('#'): continue

        list_sColumn            = sReadLine.strip('\n').split('\t')
        cDB                     = cDenovoData()
        cDB.sSampleID           = list_sColumn[0]
        cDB.sStudyName          = list_sColumn[1]
        cDB.sPubmedID           = list_sColumn[2]
        cDB.nNumProbands        = list_sColumn[3]
        cDB.nNumControls        = list_sColumn[4]
        cDB.sSequenceType       = list_sColumn[5]
        cDB.sPrimaryPhenotype   = list_sColumn[6]
        cDB.sValidation         = list_sColumn[7]
        cDB.sChr                = list_sColumn[8]
        cDB.nPosition           = list_sColumn[9]
        cDB.sVariant            = list_sColumn[10]
        cDB.sRSID               = list_sColumn[11]
        cDB.sDbsnpBuild         = list_sColumn[12]
        cDB.sAncestralAllele    = list_sColumn[13]
        cDB.n1000GenomeCount    = list_sColumn[14]
        cDB.nExacFreq           = list_sColumn[15]
        cDB.nEspAaFreq          = list_sColumn[16]
        cDB.nEspEaFreq          = list_sColumn[17]
        cDB.sTranscript         = list_sColumn[18]
        cDB.nCodingDnaSize      = list_sColumn[19]
        cDB.sGeneSym            = list_sColumn[20]
        cDB.sFunctionClass      = list_sColumn[21]
        cDB.sCDNAVariant        = list_sColumn[22]
        cDB.sProteinVariant     = list_sColumn[23]
        cDB.sExonORIntron       = list_sColumn[24]
        cDB.nPolyPhen1          = list_sColumn[25]
        cDB.nPolyPhen2          = list_sColumn[26]
        cDB.nSiftScore          = list_sColumn[27]
        cDB.nCaddScore          = list_sColumn[28]
        cDB.nLofScore           = list_sColumn[29]
        cDB.nLrtScore           = list_sColumn[30]
        list_sOutput.append(cDB)
    #loop END: sReadLine
    InFile.close()

    return list_sOutput
#def END: parse_db


def convert_to_vcf (sOutDir, list_cDB):

    sOutFile = '%s/denovodb.vcf' % sOutDir
    OutFile  = open(sOutFile, 'w')

    sHeader  = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
               % ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
    OutFile.write(sHeader)
    for cDB in list_cDB:

        sOut = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
               % (cDB.sChr, cDB.nPosition, '.', cDB.sVariant.split('>')[0], cDB.sVariant.split('>')[1],
                  '.', '.', '%s-%s' % (cDB.sSequenceType, cDB.sPrimaryPhenotype), 'GT:GQ:DP:BQ:MQ:AD:FA:VAQ:SS:FT',
                  './.:0:1:.:0:.:.:0:.:IntersectionFailure-SnpFilter')
        OutFile.write(sOut)
    #loop END: cDB
    OutFile.close()
#def END: convert_to_vcf


def convert_to_vcf_v2 (sOutDir, list_cDB):

    sOutFile = '%s/npdenovo.vcf' % sOutDir
    OutFile  = open(sOutFile, 'w')
    sHeader  = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
               % ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
    OutFile.write(sHeader)
    for cDB in list_cDB:
        sOut = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
               % (cDB.sChrID, cDB.nStartPos, '.', cDB.sRefNuc, cDB.sAltNuc,
                  '.', '.', '%s-%s' % (cDB.sSeqData, cDB.sStudy), 'GT:GQ:DP:BQ:MQ:AD:FA:VAQ:SS:FT',
                  './.:0:1:.:0:.:.:0:.:IntersectionFailure-SnpFilter')
        OutFile.write(sOut)
    #loop END: cDB
    OutFile.close()
#def END: convert_to_vcf


def annovar_vcf (sWorkDir):

    bTestRun  =  True
    sQueue    = 'optiplex.q'
    sJobName  = 'Jinman.ANNOVAR.VCF.DenovoDB'
    sLogDir   = '%s/log/%s/%s' % (sBASE_DIR, sJobName, sTIME_STAMP)
    sTempDir  = '%s/temp' % sWorkDir
    os.makedirs(sTempDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    sANNOVAR  = '/extdata6/Jinman/util/annovar/table_annovar.pl'
    sVCFFile  = '%s/denovodb.vcf'      % sWorkDir
    sAnnoFile = '%s/denovodb.anno.vcf' % sWorkDir

    # Annotate with ANNOVAR
    sScript   = '%s %s '            % (sANNOVAR, sVCFFile)
    sScript   += '/extdata6/Daekwan/util/annovar/humandb '
    sScript   += '-buildver hg19 '
    sScript   += '--vcfinput '
    sScript   += '--outfile %s '    % sAnnoFile
    sScript   += '--tempdir %s '    % sTempDir
    sScript   += '-protocol refGene,ensGene '
    sScript   += '-operation  g,g, ;'

    if bTestRun:
        print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s'
                  % (sScript, sLogDir, sQueue, sJobName))
#def END: annovar_vcf


def annovar_vcf_v2 (sWorkDir):

    bTestRun  =  False
    sQueue    = 'optiplex.q'
    sJobName  = 'Jinman.ANNOVAR.VCF.NPDenovo'
    sLogDir   = '%s/log/%s/%s' % (sBASE_DIR, sJobName, sTIME_STAMP)
    sTempDir  = '%s/temp' % sWorkDir
    os.makedirs(sTempDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    sANNOVAR  = '/extdata6/Jinman/util/annovar/table_annovar.pl'
    sVCFFile  = '%s/npdenovo.vcf'      % sWorkDir
    sAnnoFile = '%s/npdenovo.vcf'      % sWorkDir

    # Annotate with ANNOVAR
    sScript   = '%s %s '            % (sANNOVAR, sVCFFile)
    sScript   += '/extdata6/Daekwan/util/annovar/humandb '
    sScript   += '-buildver hg19 '
    sScript   += '--vcfinput '
    sScript   += '--outfile %s '    % sAnnoFile
    sScript   += '--tempdir %s '    % sTempDir
    sScript   += '-protocol refGene,ensGene '
    sScript   += '-operation  g,g, ;'

    if bTestRun:
        print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N %s'
                  % (sScript, sLogDir, sQueue, sJobName))
#def END: annovar_vcf


def parse_db_annovar (sInFile):

    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:
        if sReadLine.startswith('#'): continue
        list_sColumn       = sReadLine.strip('\n').split('\t')
        cDB                = cDenovoData()
        cDB.sChrID         = list_sColumn[0]
        cDB.nPos           = int(list_sColumn[1])
        cDB.sDBSNP_ID      = list_sColumn[2]
        cDB.sRefNuc        = list_sColumn[3]
        cDB.sAltNuc        = list_sColumn[4]
        cDB.sFilter        = list_sColumn[5]
        cDB.sQual          = list_sColumn[6]
        cDB.sInfo          = list_sColumn[7]

        # Split the Info column by ';' then key:value by '='
        dict_sInfo   = dict(sInfo.split('=') for sInfo in cDB.sInfo.split(';')[:-1] if len(sInfo.split('=')) == 2)
        cDB.sLocale  = dict_sInfo['Func.refGene'].split('\\x3b')[0]
        cDB.sGeneSym = dict_sInfo['Gene.refGene'].split('\\x3b')[0].upper()
        cDB.sGeneID  = dict_sInfo['Gene.ensGene'].split('\\x3b')[0].upper()

        list_sOutput.append(cDB)
    #loop END: sReadLine
    InFile.close()
    return list_sOutput
#def END: parse_db_annovar

## endregion

## region DBD Database

def dbd_data ():

    sWorkDir     = '%s/ref'           % sBASE_DIR
    sDBDFile     = '%s/DBD_data.txt'  % sWorkDir
    sNPDFile     = '%s/NPdenovo.txt'  % sWorkDir

    ## DBD database
    list_cDB     = parse_dbd (sDBDFile)
    list_cDB     = [cDB for cDB in list_cDB if cDB.sInheritance == 'De novo']
    print(len(list_cDB))

    ## NP-denovo database
    list_cDB    = parse_npd (sNPDFile)
    print(len(list_cDB))

    #convert_to_vcf_v2 (sWorkDir, list_cDB)
    #annovar_vcf_v2 (sWorkDir)


    sInFile       = '%s/npdenovo.vcf.hg19_multianno.vcf' % sWorkDir
    list_cDB_anno = parse_db_annovar (sInFile)

    dict_sLocale  = {}

    for cDB in list_cDB_anno:
        if cDB.sLocale not in dict_sLocale:
            dict_sLocale[cDB.sLocale] = 0
        dict_sLocale[cDB.sLocale] += 1

    for sLocale in dict_sLocale:
        print(sLocale, dict_sLocale[sLocale])


    sys.exit()

    pass
#def END: dbd_data

def parse_dbd (sInFile):
    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:
        ''' File Format
        0	Gene    	        BRAT1
        1	ID/DD
        2	Autism                      sDisease
        3	Epilepsy	        X       sDisease
        4	ADHD	                    sDisease
        5	Schizophrenia	            sDisease
        6	Bipolar Disorder	        sDisease
        7	WGS	                        sSeqData
        8	WES	X                       sSeqData
        9	Genome-wide CMA	            sSeqData
        10	Targeted CNV analysis	    sSeqData
        11	Targeted sequencing	        sSeqData
        12	Variant Type	    Frameshift
        13	Inheritance	        Bi-parental
        14	Chr
        15	Start
        16	Stop
        17	Size (bp)
        18	Genome Build
        19	Coding DNA Change	c.962_963del
        20	Protein Change	    p.Leu321Profs*81
        21	Position
        22	Reference Sequence
        23	Alternate Sequence
        24	Reference	        Saitsu et al. 2014
        25	Individual ID
        26	Additional Information
        27	PMID	            25319849
        28	Reviewers	        2
        '''
        if sReadLine.startswith('Gene'): continue
        if '"' in sReadLine: continue
        list_sColumn            = sReadLine.strip('\n').split(',')

        cDB                     = cDenovoData()
        cDB.sGeneSym            = list_sColumn[0]

        dict_sDiseases          = {1:'IntellectDisable', 2:'Autism', 3:'Epilepsy', 4:'ADHD', 5:'Schizophrenia', 6:'BipolarDisorder'}
        cDB.sDisease            = dict_sDiseases[[i+1 for i in range(6) if list_sColumn[i+1]][0]]

        dict_sSeqData           = {7:'WGS', 8:'WES', 9:'Genome-Wide_CMA', 10:'TargetedCNV', 11:'Targeted_Seq'}
        cDB.sSeqData            = dict_sSeqData[[i+7 for i in range(5) if list_sColumn[i+7]][0]]

        cDB.sVariantType        = list_sColumn[12]
        cDB.sInheritance        = list_sColumn[13]
        cDB.sChrID              = list_sColumn[14]
        cDB.nStartPos           = int(list_sColumn[15]) if list_sColumn[15] else list_sColumn[15]
        cDB.nEndPos             = int(list_sColumn[16]) if list_sColumn[16] else list_sColumn[16]
        cDB.nSize               = int(list_sColumn[17]) if list_sColumn[17] else list_sColumn[17]
        cDB.sGenomeBuild        = list_sColumn[18]
        cDB.sSubstitution       = list_sColumn[19]
        cDB.sAAchange           = list_sColumn[20]

        cDB.nVarPos             = int(list_sColumn[21]) if list_sColumn[21] else list_sColumn[21]
        cDB.sRefNuc             = list_sColumn[22]
        cDB.sAltNuc             = list_sColumn[23]
        cDB.sStudyName          = list_sColumn[24]
        cDB.sIndyID             = list_sColumn[25]
        cDB.sAddInfo            = list_sColumn[26]
        cDB.sPMID               = list_sColumn[27]
        cDB.nReviewers          = int(list_sColumn[28])
        list_sOutput.append(cDB)
    #loop END: sReadLine
    InFile.close()
    return list_sOutput
#def END: parse_db


def parse_npd (sInFile):
    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:
        ''' File Format
        0	Chr	        chr1
        1	Start	    41976742
        2	End	        41976742
        3	Ref	        G
        4	Alt	        A
        5	Gene_symbol	HIVEP3
        6	Disease	    ASD
        7	Sequencing	WGS
        8	Study	    Michaelson JJ et al. Cell 2012
        9	PMID	    23260136

        '''
        if sReadLine.startswith('Chr'): continue
        list_sColumn            = sReadLine.strip('\n').split('\t')

        cDB                     = cDenovoData()
        cDB.sChrID              = list_sColumn[0]
        cDB.nStartPos           = int(list_sColumn[1])
        cDB.nEndPos             = int(list_sColumn[2])
        cDB.sRefNuc             = list_sColumn[3]
        cDB.sAltNuc             = list_sColumn[4]
        cDB.sGeneSym            = list_sColumn[5]
        cDB.sDisease            = list_sColumn[6]
        cDB.sSeqData            = list_sColumn[7]
        cDB.sStudy              = list_sColumn[8]
        cDB.sPMID               = list_sColumn[9]

        list_sOutput.append(cDB)
    #loop END: sReadLine
    InFile.close()

    return list_sOutput
#def END: parse_db

## endregion

def main():
    print('Function Name Required')
    pass

if __name__ == '__main__':
    if len(sys.argv) == 1: main()
    else:
        function_name       = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys(): locals()[function_name](*function_parameters)
        else: sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    #if END: len(sys.argv)
#if END: __name__