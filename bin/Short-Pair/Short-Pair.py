#!/usr/bin/env python3
import os
import sys
import pickle
import re
import subprocess
import math
from optparse import OptionParser
from numpy import percentile
from numpy import std
from numpy import mean
from os.path import abspath, dirname, realpath	# Tocado JT



sequence_length_global = 0
exec_path = abspath(dirname(realpath(__file__)))	# Tocado JT

def get_first_sequence_in_fasta_file(fasta_file):
    fi = open(fasta_file, 'rt')
    sequence = None
    for line in fi:
        if line[0] == '>':
            if sequence:
                break
            sequence = ''
        else:
            sequence += line.strip()
    fi.close()
    return sequence


def step1(fastaFile, pattern):
    command = "%s/DNA2Protein 1-6 %s %s" % (exec_path, fastaFile, pattern) # exec_path, FPS
    os.system(command)
    
def step2(hmm, inputFile):
    inputFileDir = abspath(os.path.dirname(inputFile))
    inputFile = inputFile.split('/')[-1]
    outputDir = os.path.join(inputFileDir, '../out1')
    
    outputFile = os.path.join(outputDir, inputFile +'.hmmer')
    command = '%s/../hmmer/hmmsearch -E 1000 --domtblout %s %s %s/%s > /dev/null' % (exec_path, outputFile, hmm, inputFileDir, inputFile) # exec_path, FPS 
    os.system(command)    
    return outputFile
    
def ExtractHMMER(inputFile, outputFile):
    fo = open(outputFile, 'w')
    with open(inputFile, 'r') as f:
        for line in f:
            if line == '' or line[0]=='#':
                continue
            row = line.split()
            read = row[0]
            family = row[4]; family = family[:family.rfind('.')]
            if 'frame1' in inputFile or 'frame2' in inputFile or 'frame3' in inputFile:
                direction = '+'
            elif 'frame4' in inputFile or 'frame5' in inputFile or 'frame6' in inputFile:
                direction = '-' 
            fmtString = '%s %s %s %s %s %s %s %s %s\n' % (read, family, row[13], row[11], row[15], row[16], row[17], row[18], direction)
            fo.write(fmtString)
    fo.close()

def step3(inputFile):
    inputFileDir = os.path.dirname(inputFile)
    parentDir = os.path.abspath(os.path.join(inputFileDir, os.pardir))
    outputFileDir = os.path.join(parentDir, 'Out_extracted')
    if not os.path.isdir(outputFileDir):
        os.makedirs(outputFileDir)
    outputFilename = os.path.basename(inputFile) + '.extractedHmmer'
    outputFile = os.path.join(outputFileDir, outputFilename)
    ExtractHMMER(inputFile, outputFile)
    return outputFile
    
def combineHmmer(inputFolder, outputFolder, allOutput, nDomains, namePattern):
    fa = open(allOutput, 'wt')
    for i in range(1, nDomains+1):     
#    for i in xrange(201, nDomains+1):
        mapLinesList = []
        mapScoresList = []
        readPairSet = set()
        for j in range(1, 7):
            inputFile = namePattern % (i, j)
            fi = open(inputFolder+'/'+inputFile, 'rt')
            mapReadPairToLines = {}
            mapReadPairToScore = {}
            for line in fi:
                row = line.split()
                read = row[0]; family = row[1]; score = float(row[2])
                pair = (read,family)
                if pair not in mapReadPairToLines:
                    mapReadPairToLines[pair] = []
                    mapReadPairToScore[pair] = score
                mapReadPairToLines[pair].append(line)
                if score > mapReadPairToScore[pair]:
                    mapReadPairToScore[pair] = score
                readPairSet.add(pair)
            fi.close()
            mapLinesList.append(mapReadPairToLines)
            mapScoresList.append(mapReadPairToScore)
        ### start to combine
        outputFile = inputFile + "combined.hmmer"
        outputFile = os.path.join(outputFolder, outputFile)
        fo = open(outputFile, 'wt')
        for pair in sorted(readPairSet):
            maxScore = -sys.maxsize        
            for j in range(1, 7):
                mapReadPairToLines = mapLinesList[j-1]
                mapReadPairToScore = mapScoresList[j-1]
                if pair not in mapReadPairToScore:
                    continue
                score = mapReadPairToScore[pair]
                if score > maxScore:
                    maxScore = score
                    lines = mapReadPairToLines[pair]
            for line in lines:
                fo.write(line.strip() + "\n")
                fa.write(line.strip() + "\n")
        fo.close() 
    fa.close()  
    
def step4(inputFolder, namePattern, outputFile):
    mapLinesList = []
    mapScoresList = []
    readPairSet = set()
    for j in range(1, 7):
        inputFile = namePattern % (j)
        fi = open(inputFolder+'/'+inputFile, 'rt')
        mapReadPairToLines = {}
        mapReadPairToScore = {}
        for line in fi:
            row = line.split()
            read = row[0]; family = row[1]; score = float(row[2])
            pair = (read,family)
            if pair not in mapReadPairToLines:
                mapReadPairToLines[pair] = []
                mapReadPairToScore[pair] = score
            mapReadPairToLines[pair].append(line)
            if score > mapReadPairToScore[pair]:
                mapReadPairToScore[pair] = score
            readPairSet.add(pair)
        fi.close()
        mapLinesList.append(mapReadPairToLines)
        mapScoresList.append(mapReadPairToScore)
    ### start to combine   
    fo = open(outputFile, 'wt')
    for pair in sorted(readPairSet):
        maxScore = -sys.maxsize        
        for j in range(1, 7):
            mapReadPairToLines = mapLinesList[j-1]
            mapReadPairToScore = mapScoresList[j-1]
            if pair not in mapReadPairToScore:
                continue
            score = mapReadPairToScore[pair]
            if score > maxScore:
                maxScore = score
                lines = mapReadPairToLines[pair]
        for line in lines:
            fo.write(line.strip() + "\n")            
    fo.close() 
        
    
def combineFiles(inputFile1, inputFile2, outputFile):
    with open(outputFile, 'wt') as fo:
        with open(inputFile1, 'rt') as f:
            fo.write(f.read())
        with open(inputFile2, 'rt') as f:
            fo.write(f.read())        
    
def step5(inputFile1, inputFile2, outputFile):
    combineFiles(inputFile1, inputFile2, outputFile)
    
            
def step6(inputFile, outputFile):
    global sequence_length_global
    upper = math.ceil(sequence_length_global/3)
    lower = upper - 4
    f = open(inputFile, 'rt')
    nTotalLine = 0; nBoth = 0; nLeft = 0; nRight = 0; nNone = 0    
    fo1 = open(outputFile, 'wt')                
    for line in f:                
        half = line.split()
        read = half[0][:-2]; family = half[1]
        leftLength = int(half[7])-int(half[6])+1        
        rightLength = int(half[7])-int(half[6])+1   
        if leftLength >= lower and leftLength <= upper:
            nBoth += 1
            fo1.write(line)        
        nTotalLine += 1          
    f.close()
    fo1.close()
    
class HMM:
    name = None
    length = None
    map = None
    nseq = None
    pairs = None
    
    def __lt__(self, other):
        return self.name < other.name
    def __le__(self, other):
        return self.name <= other.name    
    
def readHmm(filename):
    f = open(filename)
#     name = filename[filename.rfind('/')+1:filename.find('.hmm')]    
    hmmList = []
    hmm = None
    for line in f:
        m = re.match('ACC\s+(\w+)\.\d', line)
        if m is not None:
            hmm = HMM()
            name = m.group(1)
            hmm.name = name
            i = 1     
            continue       
        m = re.match('LENG  (\d+)', line)
        if m is not None:
            length = int(m.group(1))
            hmm.length = length
            hmm.pairs = []
            hmmList.append(hmm)
            continue
        m = re.match('MAP   (\w+)', line)
        if m is not None:
            map = m.group(1)
            hmm.map = map
            continue
        m = re.match('NSEQ  (\d+)', line)
        if m is not None:
            nseq = int(m.group(1))
            hmm.nseq = nseq
            continue
        if hmm != None and hmm.length is not None:        
            row = line.split()
            if row[0] == str(i):
                six = row[21]
                pair = (i, six)
                hmm.pairs.append(pair)
                i += 1
                continue
    f.close()
    return hmmList   
    
def step7(inputFile):
    hmmYList = []; hmmNList = []
    hmmList = readHmm(inputFile)
    for hmm in hmmList:
        if hmm.map == 'yes':
            hmmYList.append(hmm)
        elif hmm.map == 'no':
            hmmYList.append(hmm)
            
    outputDir = os.path.dirname(inputFile)
    outputFile = os.path.join(outputDir, 'hmms.sav')    
    fo = open(outputFile,'wb')
    pickle.dump(hmmYList, fo)    
    fo.close()
    return outputFile, hmmList
    
class StockholmFamily:
    name = None
    sequences = None
    seq_cons = None
    def __init__(self):
        self.sequences = [] 
            
def step8(inputFile):
    families = []
    fi = open(inputFile, 'rt')
    line = fi.readline()
    inFamily = False
    while line != '':
        ### header
        if not inFamily:        
            m = re.match("#=GF AC\s+([\.\w]+)", line)
            if m != None:
                family = m.group(1)
                pos = family.rfind('.')
                family = family[:pos]
                stockholmF =  StockholmFamily()
                stockholmF.name = family
                inFamily = True
        ### end        
        if line.strip()=='//':
            inFamily = False
            families.append(stockholmF)
        ### sequences        
        if inFamily and line[0]!='#':
            row = line.split()
            sequence = (row[0], row[1])
            stockholmF.sequences.append(sequence)
        ### consensus sequence
        if inFamily:
            m = re.match("#=GC seq_cons\s+(.+)", line)
            if m != None:
                stockholmF.seq_cons = m.group(1)
        line = fi.readline()
        
    fi.close()
    
    dumpFile = inputFile+'.sav'
    f = open(dumpFile, 'wb')
    pickle.dump(families, f)
    f.close()
    return dumpFile

def LoadReadPairs(filename):
    f = open(filename, 'rt')
    posReadSet = set(); negReadSet = set()
    mapPosReadToAlignments = {}; mapNegReadToAlignments = {} 
    familySet = set()
    for line in f:
        row = line.split()
        read = row[0][:-2]
        dir = row[0][-2:]
        familySet.add(row[1])
        if dir == '.1':
            posReadSet.add(read)
            if read not in mapPosReadToAlignments:
                mapPosReadToAlignments[read] = []
            mapPosReadToAlignments[read].append(row)
        else:
            negReadSet.add(read)
            if read not in mapNegReadToAlignments:
                mapNegReadToAlignments[read] = []
            mapNegReadToAlignments[read].append(row)
        
    f.close()
    return posReadSet & negReadSet, mapPosReadToAlignments, mapNegReadToAlignments, posReadSet, negReadSet, familySet  
    
def CalculateFragmentLengths2(filename, hmmSav, seedSav):
    ### load the hmm data
    inputFile = hmmSav
    fi = open(inputFile, 'rb')
    hmmYList = pickle.load(fi)
    fi.close()
    stockholmF = StockholmFamily()
    inputFile = seedSav
    fi = open(inputFile, 'rb')
    familyList = pickle.load(fi)
    fi.close()
    
    readPairSet, mapPosReadToAlignments, mapNegReadToAlignments, posReadSet, negReadSet, familySet = LoadReadPairs(filename);
    lengthList = []
    nNoPair = 0; nHasPair = 0; nHasMore = 0
    baseDirectory = os.path.dirname(filename)
        
    for read in readPairSet:
        alignments1 = mapPosReadToAlignments[read];
        alignments2 = mapNegReadToAlignments[read];
        if read == 'gnl|SRA|SRR360147.2138780':            
            pass
        curLength = -1; 
        record = True;
        hasPair = False
        hasMore = False
        nAlignmentCount = 0
        hasAlignment = len(alignments1) > 0 and len(alignments2) > 0
        fmtString = '%s\n' % (read)
    
        for alignment1 in alignments1:
            for alignment2 in alignments2:
                if alignment2[1] != alignment1[1]:  continue
                if not hasMore:
                    for alignment11 in alignments1:
                        if alignment11 == alignment1:   continue
                        if alignment11[1] == alignment1[1]:
                            hasMore = True
    
                            break
                if not hasMore:
                    for alignment22 in alignments2:
                        if alignment22 == alignment2:   continue
                        if alignment22[1] == alignment1[1]:
                            break 
                hasPair = True
                nAlignmentCount += 2
                
                ### new way to calculate fragment length
                a1 = int(alignment1[5]); b1 = int(alignment1[4])
                a2 = int(alignment2[5]); b2 = int(alignment2[4])
                for hmm in hmmYList:
                    if hmm.name == alignment1[1]:
                        break
                stockholmF = StockholmFamily()
                for stockholmF in familyList:
                    if stockholmF.name == alignment1[1]:
                        break 
                try:            
                    aa1 = int(hmm.pairs[a1-1][1])
                    bb1 = int(hmm.pairs[b1-1][1])
                    aa2 = int(hmm.pairs[a2-1][1])
                    bb2 = int(hmm.pairs[b2-1][1])
                except:                    
                    continue
                
                def calculateNewPos(a, sequence):
                    sequence = sequence[1]
                    x1 = sequence[:a-1]
                    nDots = len(x1)-len(x1.replace('.',''))
                    na1 = a-nDots  
                    return na1
                lengthList2 = []                                  
                for sequence in stockholmF.sequences:
                    na1 = calculateNewPos(aa1, sequence)
                    nb1 = calculateNewPos(bb1, sequence)
                    na2 = calculateNewPos(aa2, sequence)
                    nb2 = calculateNewPos(bb2, sequence)
                    length = max(na1, na2) - min(nb1, nb2) + 1
                    lengthList2.append(length)
                length = mean(lengthList2)
                ### end new way to calcualte fragment length
                if read == 'gnl|SRA|SRR360147.2138780':
                    continue                
                if curLength == -1: curLength = length;
                
                if curLength != length:
                    record = False
                    break
            if not record:
                break
        if hasPair and record:
            if read == 'gnl|SRA|SRR360147.2138780':                
                continue
                
            lengthList.append(curLength)
    
        if hasAlignment and not hasPair:
            nNoPair += 1
        if hasPair:
            nHasPair += 1
        if hasMore:
            nHasMore += 1    
    
    outputFile = '%s/fragment_length_glocal1.txt' % (baseDirectory)
    f = open(outputFile, 'wt')
    for length in lengthList:
        f.write('%d\n' % length)
    f.close()    
    return outputFile
    
def step9(filename, hmmSav, seedSav):
    return CalculateFragmentLengths2(filename, hmmSav, seedSav)
    
def step10(inputFile):
    baseDirectory = os.path.dirname(inputFile)
    
    lengthList = []
    f = open(inputFile, 'rt')
    for length in f:
        length = int(length)
        lengthList.append(length)
    f.close()
    
    mapLengthToNumber = {}
    for length in lengthList:
        if length not in mapLengthToNumber:
            mapLengthToNumber[length] = 0
        mapLengthToNumber[length] += 1
    outputFile = '%s/fragment_length_dist_glocal1.txt' % (baseDirectory)
    f = open(outputFile, 'wt')
    for length in sorted(mapLengthToNumber.keys()):
        fmtString = '%d,%d\n' % (length, mapLengthToNumber[length])
        f.write(fmtString)
    f.close()  
    return outputFile   

def step11(inputFile, outputFile):
    posReadFamilySet = set(); negReadFamilySet = set()
    mapPosReadFamilyToLine = {}; mapNegReadFamilyToLine = {}
    with open(inputFile, 'rt') as fi:
        for line in fi:
            row = line.split()
            name = row[0][:-2]; dir = row[0][-1]; family = row[1]
            if dir=='1':
                posReadFamilySet.add((name, family))
                mapPosReadFamilyToLine[(name,family)] = line
            elif dir=='2':
                negReadFamilySet.add((name, family))
                mapNegReadFamilyToLine[(name,family)] = line
    with open(outputFile, 'wt') as fo:
        for (name, family) in sorted(posReadFamilySet):
            fmtString = mapPosReadFamilyToLine[(name,family)].rstrip()
            if (name,family) in negReadFamilySet:
                fmtString = '%s,%s' % (fmtString, mapNegReadFamilyToLine[(name,family)].rstrip())
            fmtString += '\n'
            fo.write(fmtString)
        for (name, family) in sorted(negReadFamilySet-posReadFamilySet):
            fmtString = mapNegReadFamilyToLine[(name,family)].rstrip()
            fmtString += '\n'
            fo.write(fmtString)
            
def step12(hmmFile, hmmList, outputFolder):
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
    for hmm in sorted(hmmList, key=lambda x:x.name):
        family = hmm.name
        command = "%s/get_hmm.sh %s %s > %s/%s.hmm" % (exec_path, family, hmmFile, outputFolder,family) # exec_path, FPS 
        os.system(command)    

def runHmmer(inputFile1, inputFile2, inputFile3, hmmFolder, outputFileInput, fastaPath, faaPath, filterFile):
    ### create the directories if they do not exist
    if not os.path.exists(fastaPath):
        os.mkdir(fastaPath)
    if not os.path.exists(faaPath):
        os.mkdir(faaPath)
    
    mapReadLeftToLine = {}

    readSet = set()
    fi = open(inputFile1, 'rt')
    fo = open(filterFile, 'wt')
    for line in fi:
        if line.find(",") >= 0:
            fo.write(", ".join(line.split(",")))
            continue
        row = line.split()
        pos = row[0].rfind(".")
        read = row[0][:pos]
        family = row[1]
        key = "%s_%s" % (row[0], family)
        readSet.add((row[0], family, key))
        if key not in mapReadLeftToLine:
            mapReadLeftToLine[key] = []
        mapReadLeftToLine[key].append(line.strip())                
    fi.close()
    fo.close()
    fo = None    
    
    readNames = set()
    for (read, family, key) in readSet:
        readNames.add(read.replace('|','_'))
    write = False   
    fi = open(inputFile2, 'rt')
    nLine = 0
    for line in fi:
        if line[0]=='>':
            read = line[1:].strip();
            read = read.replace("|","_")
            if read[:-2]+".2" in readNames:
                outputFile = r'%s/%s.fasta' % (fastaPath, read)
                fo = open(outputFile, 'wt')
                fo.write(line)
                write = True
        else:
            if write:               
                fo.write(line)
                fo.close()
                write = False
        nLine += 1
    fi.close()           
    write = False    
    
    fi = open(inputFile3, 'rt')
    for line in fi:
        if line[0]=='>':
            read = line[1:].strip();
            read = read.replace("|","_")            
            if read[:-2]+".1" in readNames:
                outputFile = r'%s/%s.fasta' % (fastaPath, read)
                fo = open(outputFile, 'wt')
                fo.write(line)
                write = True
        else:
            if write:
                fo.write(line)
                fo.close()
                write = False
    fi.close()
    

    foA = open(outputFileInput, 'wt')
    for (read, family, key) in sorted(readSet):
        tmpRead =read.replace("|","_")
        if tmpRead[-2:] == '.1':
            tmpRead = tmpRead[:-2] + '.2'
        elif tmpRead[-2:] == '.2':
            tmpRead = tmpRead[:-2] + '.1'
        fastaFile = r'%s/%s.fasta' % (fastaPath, tmpRead)
        faaFile = r'%s/%s.faa' % (faaPath, tmpRead)
        command = "%s/DNA2Protein 1-6 %s %s" % (exec_path, fastaFile, faaFile) # exec_path, FPS
        os.system(command)
        outputFile = r'%s.allframe' % (faaFile)           
        fo = open(outputFile, 'wt')
        maxScore = -sys.maxsize
        maxString = ''
        for i in range(1,7):
            command = "%s/hmmer3_pipeline_missing_end.sh" % exec_path # exec_path, FPS
            args = [command, "%s/%s.hmm" % (hmmFolder, family), "%s.frame%d" % (faaFile, i), "%d" % i, "10"]
            my_env = os.environ # Add SQM/bin/hmmer to the PATH env. variable for subprocess
            my_env['PATH'] = '%s/../hmmer:' % exec_path + os.environ['PATH']
            outString = subprocess.Popen(args, env=my_env, stdout=subprocess.PIPE).communicate()[0].decode()
            fo.write(outString)
            if outString:
                row = outString.split()
                score = float(row[2])
                if score > maxScore:
                    maxScore = score
                    maxString = outString
        fo.close()
        for line in mapReadLeftToLine[key]:
            fmtString = "%s, %s\n" % (line.strip(), maxString.strip())
            foA.write(fmtString)
        
    foA.close()

def step13(inputFile1, inputFile2, inputFile3, hmmFolder, outputFileInput, fastaPath, faaPath, filterFile):
    runHmmer(inputFile1, inputFile2, inputFile3, hmmFolder, outputFileInput, fastaPath, faaPath, filterFile)
    
def removeOneEnd(inputFile, outputFile):
    fi = open(inputFile, 'rt')    
    fo = open(outputFile, 'wt')
    nBoth = 0; nTotal = 0
    for line in fi:
        nTotal += 1        
        if line.rstrip()[-1]==',':
            continue
        else:
            row = line.rstrip().rstrip(',').split(',')
            if len(row)>1:
                fo.write(line)
                nBoth += 1                
    fo.close()
    fi.close()
    
def step14(inputFile, outputFile):
    removeOneEnd(inputFile, outputFile)
    
def step15(inputFile1, inputFile2, outputFile):
    combineFiles(inputFile1, inputFile2, outputFile)
    
def CalculateLikelihood2(pair, hmmYList, familyList):
    global mapLengthToFrequency
    
    row = pair[0].split()
    alignment1 = pair[0].split()
    score1 = float(row[2])       
    row = pair[1].split()
    alignment2 = pair[1].split()
    score2 = float(row[2])    
    ### new way to calculate fragment length
    a1 = int(alignment1[5]); b1 = int(alignment1[4])
    a2 = int(alignment2[5]); b2 = int(alignment2[4])
    for hmm in hmmYList:
        if hmm.name == alignment1[1]:
            break
    stockholmF = StockholmFamily()
    for stockholmF in familyList:
        if stockholmF.name == alignment1[1]:
            break 
    try:            
        aa1 = int(hmm.pairs[a1-1][1])
        bb1 = int(hmm.pairs[b1-1][1])
        aa2 = int(hmm.pairs[a2-1][1])
        bb2 = int(hmm.pairs[b2-1][1])
    except:
        pass
    
    def calculateNewPos(a, sequence):
        sequence = sequence[1]
        x1 = sequence[:a-1]
        nDots = len(x1)-len(x1.replace('.',''))
        na1 = a-nDots  
        return na1
    lengthList2 = []                                  
    for sequence in stockholmF.sequences:
        na1 = calculateNewPos(aa1, sequence)
        nb1 = calculateNewPos(bb1, sequence)
        na2 = calculateNewPos(aa2, sequence)
        nb2 = calculateNewPos(bb2, sequence)
        length = max(na1, na2) - min(nb1, nb2) + 1
        lengthList2.append(length)
    fragmentLength = int(math.floor(mean(lengthList2)))
    if fragmentLength in mapLengthToFrequency:
        pd = mapLengthToFrequency[fragmentLength]
    else:
        pd = 0.0 

    read = row[0][:-2]
    try:
        y = score1 + score2 + math.log(pd, 2)
    except Exception:
        return 0 , read, row[1]
    return y , read, row[1]    
    
def PairLikelihoodV2(hmmSav, seedSav, fragmentLengthDist, inputFile4, threshold, outputFile2):
    ### load the hmm data
    inputFile = hmmSav
    fi = open(inputFile, 'rb')
    hmmYList = pickle.load(fi)
    fi.close()
    
    inputFile = seedSav
    fi = open(inputFile, 'rb')
    familyList = pickle.load(fi)
    fi.close()
    
    global mapLengthToFrequency
    mapLengthToFrequency = {}
    
    inputFile = fragmentLengthDist
    
    total = 0
    fi = open(inputFile, 'rt')
    for line in fi:
        row = line.split(",")
        length = int(row[0])
        frequency = int(row[1])
        mapLengthToFrequency[length] = frequency
        total += frequency
    fi.close()
    for key in list(mapLengthToFrequency.keys()):
        mapLengthToFrequency[key] = mapLengthToFrequency[key] / float(total)    
    
    mapReadToLine = {}    
    mapReadToFamilySet = {}
    mapReadFamilyToLikelihood = {}
    inputFile = inputFile4
#     outputFile = outputFile1
#     fo = open(outputFile, 'wt')
    fi = open(inputFile, 'rt')
    count = 0
    for line in fi:
        pair = line.split(",")
        likelihood, read, family = CalculateLikelihood2(pair, hmmYList, familyList)
        fmtString = "%s\t%g\n" % (line.strip(), likelihood)
#         fo.write(fmtString)
        if read not in mapReadToLine:
            mapReadToLine[read] = fmtString
        if likelihood > float(fmtString.split()[-1]):
            mapReadToLine[read] = fmtString 
        if read not in mapReadToFamilySet:
            mapReadToFamilySet[read] = set()
        mapReadToFamilySet[read].add(family)
        if read not in mapReadFamilyToLikelihood:
            mapReadFamilyToLikelihood[read] = {}
        if family not in mapReadFamilyToLikelihood[read]:
            mapReadFamilyToLikelihood[read][family] = fmtString
        if likelihood > float(mapReadFamilyToLikelihood[read][family].split()[-1]):
            mapReadFamilyToLikelihood[read][family] = fmtString
        count += 1      
    fi.close()
#     fo.close()    
    
    outputFile = outputFile2
    fo = open(outputFile, 'wt')
    for read in sorted(mapReadFamilyToLikelihood.keys()):
        maxLikelihood = -1
        likelihoodArray = []
        allNeg = True
        mapLikelihoodToFamily = {}
        for family in mapReadFamilyToLikelihood[read]:
            a = float(mapReadFamilyToLikelihood[read][family].split()[-1])
            mapLikelihoodToFamily[a] = family
            if a >= 0: allNeg = False  
            likelihoodArray.append(a)
        if allNeg:
            newLikelihoodArray = []
            for a in likelihoodArray:
                a = abs(a)
                newLikelihoodArray.append(a)
                mapLikelihoodToFamily[a] = mapLikelihoodToFamily[-a]
                if a > maxLikelihood:
                    maxLikelihood = a
            for a in newLikelihoodArray:
                if a >= maxLikelihood * threshold:
                    family = mapLikelihoodToFamily[a]
                    fmtString = mapReadFamilyToLikelihood[read][family]
                    fmtString = ' '.join(fmtString.split()[:-1]) + '\n'
                    fo.write(fmtString)
        else:
            for a in likelihoodArray:
                if a > maxLikelihood:
                    maxLikelihood = a
            for a in likelihoodArray:
                if a >= maxLikelihood * threshold:
                    family = mapLikelihoodToFamily[a]
                    fmtString = mapReadFamilyToLikelihood[read][family]
                    fmtString = ' '.join(fmtString.split()[:-1]) + '\n'
                    fo.write(fmtString)       
    fo.close()     
    
    mapNumberToFrequency = {}
    for read in list(mapReadToFamilySet.keys()):
        number = len(mapReadToFamilySet[read])
        if number not in mapNumberToFrequency:
            mapNumberToFrequency[number] = 0
        mapNumberToFrequency[number] += 1
#     outputFile = outputFile3
#     fo = open(outputFile, 'wt')
#     for number in sorted(mapNumberToFrequency.keys()):
#         fo.write("%d\t%d\n" % (number, mapNumberToFrequency[number]))
#     fo.close() 
                 
def step16(hmmSav, seedSav, fragmentLengthDist, inputFile4, threshold, outputFile2):
    PairLikelihoodV2(hmmSav, seedSav, fragmentLengthDist, inputFile4, threshold, outputFile2)                 

def part1(fastaName, fastaFile1, fastaFile2, pattern1, pattern2, hmm, step4OutputFile1, step4OutputFile2, step5OutputFile):
    global sequence_length_global
    sequence_length_global = len(get_first_sequence_in_fasta_file(fastaFile1))
    step1(fastaFile1, pattern1)
    step1(fastaFile2, pattern2)
    step3InputList = []
    step3OutputList = []
    for i in range(1,7):
        inputFile = pattern1 + '.frame%d' % i
        step2Output = step2(hmm, inputFile)
        step3InputList.append(step2Output)
    for i in range(1,7):
        inputFile = pattern2 + '.frame%d' % i
        step2Output = step2(hmm, inputFile)
        step3InputList.append(step2Output)        
    for step3Input in step3InputList:        
        step3OutputList.append(step3(step3Input))
    step3OutputFolder = os.path.dirname(step3OutputList[0])
    namePattern = fastaName + '.1.frame%d.hmmer.extractedHmmer'    
    step4(step3OutputFolder, namePattern, step4OutputFile1)
    namePattern = fastaName + '.2.frame%d.hmmer.extractedHmmer'
    step4(step3OutputFolder, namePattern, step4OutputFile2)        
    step5(step4OutputFile1, step4OutputFile2, step5OutputFile)
    
def part2(inputFile, hmmFile, seedFile):
    outputFile = inputFile + '.out'
    step6(inputFile, outputFile)
    hmmSav, hmmList = step7(hmmFile)
    seedSav = step8(seedFile)
    step9OutputFile = step9(outputFile, hmmSav, seedSav)
    step10OutputFile = step10(step9OutputFile)
    return hmmSav, hmmList, seedSav, step10OutputFile
    
def part3(fastaFile1, fastaFile2, step5Output, hmmFile, hmmSav, hmmList, seedSav, hmmFolder, fragmentLengthDist, threshold, outputFile2):    
    step11Output = step5Output+'1'
    step11(step5Output, step11Output)
    step12(hmmFile, hmmList, hmmFolder)
    step13OutputFile1 = step11Output+'_1'; step13OutputFile2 = step11Output+'_2'
    workingDir = os.path.dirname(step11Output)
    fastaPath = '%s/fastaSP' % workingDir
    faaPath = '%s/faaSP' % workingDir
    step13(step11Output, fastaFile1, fastaFile2, hmmFolder, step13OutputFile1, fastaPath, faaPath, step13OutputFile2)
    step14Output = step13OutputFile1 + '.remove'
    step14(step13OutputFile1, step14Output)
    step15Output = step11Output+'_1_2'
    step15(step14Output, step13OutputFile2, step15Output)
    step16(hmmSav, seedSav, fragmentLengthDist, step15Output, threshold, outputFile2)
    
def control(fastaFile1, fastaFile2, hmmFile, seedFile, threshold, outputFile):
    baseFile = os.path.basename(fastaFile1)
    fastaName = baseFile.replace('.1.fasta','')
    baseDir = abspath(os.path.dirname(fastaFile1))
    patternDir = '%s/Protein' % baseDir
    
    if not os.path.isdir(patternDir):
        os.makedirs(patternDir)
    pattern1 = '%s/%s' % (patternDir, fastaName+'.1')
    pattern2 = '%s/%s' % (patternDir, fastaName+'.2')
    
    
    step4OutputFile1 = fastaFile1+'.alldomains.allframe.hmmer'
    step4OutputFile2 = fastaFile2+'.alldomains.allframe.hmmer'
    step5OutputFile = os.path.join(baseDir, '%s.both.alldomains.allframe.hmmer' % fastaName)
    part1(fastaName, fastaFile1, fastaFile2, pattern1, pattern2, hmmFile, step4OutputFile1, step4OutputFile2, step5OutputFile)
    hmmSav, hmmList, seedSav, step10OutputFile = part2(step5OutputFile, hmmFile, seedFile)
    hmmFolder = os.path.join(baseDir, 'HMMs')
    part3(fastaFile1, fastaFile2, step5OutputFile, hmmFile, hmmSav, hmmList, seedSav, hmmFolder, step10OutputFile, threshold, outputFile)
    
if __name__=='__main__':
    
    parser = OptionParser()
    parser.add_option("-m", "--hmm", dest="hmmFile",
                      help="HMMER3 HMM file", metavar="FILE")
    parser.add_option("-s", "--seed", dest="seedFile", action="store",
                      help="seed file corresponding to HMM file")
    parser.add_option("-x", "--fasta1", dest="fastaFile1", action="store",
                      help="fasta file1")
    parser.add_option("-y", "--fasta2", dest="fastaFile2", action="store",
                  help="fasta file2")
    parser.add_option("-t", "--threshold", dest="threshold", action="store",
                      type="float", default=0.4, help="threashold default 0.4")
    parser.add_option("-o", "--output", dest="outputFile", action="store",
                      type="string", help="output file ")
    
     
    (options, args) = parser.parse_args()    
    
    control(options.fastaFile1, options.fastaFile2, options.hmmFile, options.seedFile, options.threshold, options.outputFile)

            
