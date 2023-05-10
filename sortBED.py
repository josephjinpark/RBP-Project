#!/extdata6/Jinman/opt/python3/bin/python3

import os, sys


def main():
    sBedTools = '/home/lab/bin/bedtools2/bin'
    sWorkDir  = '/extdata6/Jinman/07_rbp'
    sBedDir   = '%s/eclipdata/bed'           % sWorkDir
    sOutDir   = '%s/eclipdata/sorted_bed'    % sWorkDir
    sLogDir   = '%s/log/SortBed'             % sWorkDir

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    sQueue    = 'optiplex.q'
    bTestRun  = False

    list_sBedFiles = os.listdir(sBedDir)
    for sBedFile in list_sBedFiles:

        sInFile  = '%s/%s' % (sBedDir, sBedFile)
        sOutFile = '%s/%s' % (sOutDir, sBedFile)


        sScript  = 'sort -k1,1 -k2,2n -k3,3n -i %s -o %s' % (sInFile, sOutFile)

        sScript  = 'echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.SortBed.%s' \
                   % (sScript, sLogDir, sQueue, sBedFile)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    #loop END: sBedFile
#def END: main

main()


