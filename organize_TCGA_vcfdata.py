#!/extdata6/Jinman/opt/python3/bin/python3

## For parsing manifest from TCGA and creating new patient IDs with softlinks to original vcf files

import os, sys, json, time
sTIME_STAMP    = '%s' % (time.ctime().replace(' ', '-').replace(':', '_'))
sBASE_DIR      = '/extdata6/Jinman/07_rbp'
sSTUDY         = 'TCGA_ESCA'
sVarCaller     = 'VCFs'
sSeqData       = 'WGS'
sTCGA_DIR      = '%s/vardata/%s/%s' % (sBASE_DIR, sSTUDY, sSeqData)
sFILEID_DIR    = '%s/00_fileids'    % sBASE_DIR
## region Class TCGA_Summary
class TCGA_Summary: pass
## endregion

def main():

    bTestRun        = False
    sQueue          = 'optiplex.q'

    sSummaryFile    = '%s/%s_%s_%s.txt'        % (sTCGA_DIR, sSTUDY, sSeqData, sVarCaller)
    list_sData      = load_summary_data (sSummaryFile)

    sJsonFile       = '%s/%s_%s_%s.json'       % (sTCGA_DIR, sSTUDY, sSeqData, sVarCaller)
    dict_sBarcode   = load_json_data (sJsonFile)
    sOutDir         = '%s/VCFs'                % sTCGA_DIR

    sJobName        = 'Jinman.OrganizeTCGA.%s' % sSTUDY
    sLogDir         = '%s/log/%s/%s'           % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    OutFile         = open('%s/%s_patidkey.txt' % (sTCGA_DIR, sSTUDY), 'w')      # Make Patient ID key file
    OutFile2        = open('%s/%s_%s.txt'       % (sFILEID_DIR, sSTUDY, sSeqData), 'w')# Make File ID file
    for i, cSum in enumerate(list_sData):


        i = i+1

        sOriginalVCF = '%s/RawVCFs/%s/%s'  % (sTCGA_DIR, cSum.sBarcode, cSum.sFilename)
        sInFile      = '%s/RawVCFs/%s/%s'  % (sTCGA_DIR, cSum.sBarcode, cSum.sFilename.replace('.gz', ''))
        sPatID       = '%s_%s_%03d'        % (sSTUDY, sVarCaller, i)
        sOutFile     = '%s/%s.gz'          % (sOutDir, sPatID)

        sOut         = '%s\t%s\t%s\n'      % (sPatID, cSum.sFilename, dict_sBarcode[cSum.sFilename])
        OutFile.write(sOut)
        OutFile2.write('%s\n' % sPatID)
        if sOriginalVCF.endswith('.gz'):

            if os.path.isfile(sOriginalVCF):
                sCmd         = 'gunzip %s; '   % sOriginalVCF
            else:
                sCmd         = ''

        else:
            sCmd         = ''
        sCmd        += 'bgzip -c %s > %s;' % (sInFile, sOutFile)
        sCmd        += 'tabix -p vcf %s;'  % sOutFile

        sScript      = 'echo "%s" | qsub -V -q %s -cwd -j y -o %s -N %s.%s.%s ' \
                       % (sCmd, sQueue, sLogDir, sJobName, sSeqData, sPatID)
        if bTestRun: print(sScript)
        else:        os.system(sScript)
    #loop END: i, cSum
    OutFile.close()
    OutFile2.close()

#def END: main

def load_json_data (sInFile):

    InFile         = open(sInFile, 'r')
    list_sJSONData = json.load(InFile)
    InFile.close()

    dict_sTissue   = {'01A':'T','01B':'T', # Primary Solid Tumor

                      '03A':'T','03B':'T', # Primary blood-derived Cancer
                      '10A':'N','10B':'N', # Blood-derived Normal
                      '11A':'N','11B':'N'} # Solid-tissue Normal

    dict_sBarcode  = {}
    dict_sTest     = {}
    for dict_sData in list_sJSONData:

        '''
        #print(dict_sData['data_category'])
        for e in dict_sData:
            print(e, dict_sData[e])
        print('******************************')

        for e in dict_sData['cases'][0]:
            print(e, dict_sData['cases'][0][e])
        print('******************************')

        for e in dict_sData['cases'][0]['samples'][0]:
            print(e, dict_sData['cases'][0]['samples'][0][e])
        '''

        sKey     = dict_sData['file_name']
        sValue   = dict_sData['associated_entities'][0]['entity_submitter_id'].split('-')[2]
        sBarCode = dict_sData['associated_entities'][0]['entity_submitter_id']
        #print(dict_sData['associated_entities'][0]['entity_submitter_id'], sKey)

        try: sTissue  = dict_sData['associated_entities'][0]['entity_submitter_id'].split('-')[3]
        except IndexError:sTissue = 'Not'
        #if dict_sTissue[sTissue] != 'T': continue


        if sTissue not in dict_sTest:
            dict_sTest[sTissue] = 0
        dict_sTest[sTissue] += 1

        if sKey not in dict_sBarcode:
            dict_sBarcode[sKey] = ''
        dict_sBarcode[sKey] = sValue
    #loop END: dict_sData


    #VS-Check
    if not dict_sBarcode:
        sys.exit('Invalid Dictionary : load_json_data : dict_sBarcode size =%d' % len(dict_sBarcode))

    return dict_sBarcode
#def END: load_json_data


def load_summary_data (sInFile):

    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for i,sReadLine in enumerate(InFile):

        if sReadLine.startswith('id'):
            #list_sHeaders        = sReadLine.strip('\n').split('\t')
            continue
        list_sColumn             = sReadLine.strip('\n').split('\t')

        #for sHeader,sExample in zip(list_sHeaders, list_sColumn):
        #    print('%s\t%s' % (sHeader, sExample))

        cSum                     = TCGA_Summary()
        cSum.sBarcode	         = list_sColumn[0]
        cSum.sFilename		     = list_sColumn[1]
        cSum.sChecksum		     = list_sColumn[2]
        cSum.sFilesize		     = list_sColumn[3]
        cSum.sState	             = list_sColumn[4]

        if 'indel' in cSum.sFilename: continue

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