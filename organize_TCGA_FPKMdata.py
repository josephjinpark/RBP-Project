#!/extdata6/Jinman/opt/python3/bin/python3

## For determining the Tumor-Normal Pair Data availability on TCGA for WGS, WES, mRNA-seq, and sRNA-seq


import os, sys, json

sBASE_DIR      = '/extdata6/Jinman/07_rbp'
sSTUDY         = 'TCGA_LAML'
sSeqData       = 'WXS'
sTCGA_DIR      = '%s/%s/%s'      % (sBASE_DIR, sSTUDY, sSeqData)


## region Class TCGA_Summary
class TCGA_Summary: pass
## endregion

def main():

    sJsonFile     = '%s/%s_FPKMfiles.json'      % (sTCGA_DIR, sSTUDY)
    sKeyFile      = '%s/%s_patidkey.txt'        % (sTCGA_DIR, sSTUDY)
    dict_sBarcode = load_json_data(sJsonFile)
    dict_sKey     = load_keyfile (sKeyFile)

    sOutFile      = '%s/%s_rnaseq-UQ.txt'       % (sTCGA_DIR, sSTUDY)
    OutFile       = open(sOutFile, 'w')
    for sBarcode in dict_sBarcode:

        if sSTUDY == 'TCGA_LIHC':

            if len(dict_sBarcode[sBarcode]) == 1: continue
            sOut = '%s\t%s\t%s\t%s\n' % \
                   (sBarcode, dict_sBarcode[sBarcode]['T'], dict_sBarcode[sBarcode]['N'], dict_sKey[sBarcode])
        else:
            try:
                sOut = '%s\t%s\t%s\t%s\n' % \
                   (sBarcode, 'None', dict_sBarcode[sBarcode]['N'], dict_sKey[sBarcode])
            except KeyError: continue

            print(sOut[:-1])

            OutFile.write(sOut)
    #loop END: sBarcode
    OutFile.close()


def load_keyfile (sInFile):

    InFile       = open(sInFile, 'r')
    list_sOutput = [sReadLine.strip('\n').split('\t') for sReadLine in InFile]
    InFile.close()
    dict_sOutput = {}

    for sPatID, sFileName, sBarcode in list_sOutput:
        if sBarcode not in dict_sOutput:
            dict_sOutput[sBarcode] = ''
        dict_sOutput[sBarcode] = sPatID
    #loop END: sPatID, sFileName, sBarcode

    #VS-Check
    if not dict_sOutput:
        sys.exit('Invalid Dictionary : load_keyfile : dict_sOutput size =%d' % len(dict_sOutput))

    return dict_sOutput
#def END: load_keyfile


def load_json_data (sInFile):
    InFile         = open(sInFile, 'r')
    list_sJSONData = json.load(InFile)
    InFile.close()

    dict_sBarcode = {}
    dict_sTissue  = {'01A':'T','01B':'T','11A':'N', '03A':'N', '03B':'N'}

    dict_sTest    = {}

    for dict_sData in list_sJSONData:

        if not dict_sData['file_name'].endswith('FPKM-UQ.txt.gz'): continue

        sKey      = dict_sData['associated_entities'][0]['entity_submitter_id'].split('-')[2]
        sTissue   = dict_sTissue[dict_sData['associated_entities'][0]['entity_submitter_id'].split('-')[3]]
        sFileName = dict_sData['file_name']

        sTestKey  = dict_sData['associated_entities'][0]['entity_submitter_id'].split('-')[3]

        if sTestKey not in dict_sTest:
            dict_sTest[sTestKey] = 0
        dict_sTest[sTestKey] += 1

        if sKey not in dict_sBarcode:
            dict_sBarcode[sKey] = {}

        if sTissue not in dict_sBarcode[sKey]:
            dict_sBarcode[sKey][sTissue] = ''
        dict_sBarcode[sKey][sTissue] = sFileName
    #VS-Check
    if not dict_sBarcode:
        sys.exit('Invalid Dictionary : load_json_data : dict_sBarcode size =%d' % len(dict_sBarcode))

    return dict_sBarcode
#def END: load_json_data


def load_summary_data (sInFile):

    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for i,sReadLine in enumerate(InFile):

        if sReadLine.startswith('study'):
            #list_sHeaders        = sReadLine.strip('\n').split('\t')
            continue
        list_sColumn             = sReadLine.strip('\n').split('\t')

        #for sHeader,sExample in zip(list_sHeaders, list_sColumn):
        #    print('%s\t%s' % (sHeader, sExample))

        cSum                     = TCGA_Summary()

        cSum.sStudy	             = list_sColumn[0]
        cSum.sBarCode	         = list_sColumn[1]
        cSum.sDisease	         = list_sColumn[2]
        cSum.sDiseaseName	     = list_sColumn[3]
        cSum.sSampleType	     = list_sColumn[4]
        cSum.sSampleTypeName	 = list_sColumn[5]
        cSum.sAnalyteType        = list_sColumn[6]
        cSum.sLibraryType	     = list_sColumn[7]
        cSum.sCenter		     = list_sColumn[8]
        cSum.sCenterName         = list_sColumn[9]
        cSum.sPlatform		     = list_sColumn[10]
        cSum.sPlatformName       = list_sColumn[11]
        cSum.sAssembly		     = list_sColumn[12]
        cSum.sFileName		     = list_sColumn[13]
        cSum.sFileSize		     = list_sColumn[14]
        cSum.sChecksum		     = list_sColumn[15]
        cSum.sAnalysisID		 = list_sColumn[16]
        cSum.sAliquotID		     = list_sColumn[17]
        cSum.sPatientID	         = list_sColumn[18]
        cSum.sSampleID		     = list_sColumn[19]
        cSum.sTssID		         = list_sColumn[20]
        cSum.sSampleAccession	 = list_sColumn[21]
        cSum.sPublished		     = list_sColumn[22]
        cSum.sUploaded		     = list_sColumn[23]
        cSum.sModified		     = list_sColumn[24]
        cSum.sState	             = list_sColumn[25]

        list_sOutput.append(cSum)
    #loop END: sReadLine

    #V-S Check:
    if not list_sOutput:
        sys.exit('Invalid List : load_summary_data : list_sOutput size= %d' % len(list_sOutput))

    return list_sOutput
#def END: load_summary_data


def characterize_TN_pairs (dict_sData):

    #list of sPatientIDs
    list_T_WGS = [cSum.sPatientID for cSum in dict_sData['Tumor']['WGS']]

    print(list_T_WGS)
#def END: characterize_TN_pairs


main()