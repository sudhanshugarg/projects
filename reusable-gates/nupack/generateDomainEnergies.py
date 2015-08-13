import re, string, sys
import subprocess
import datetime
#from bs4 import BeautifulSoup
import requests

def runMfe(which):
   subprocess.call(['mfe', '-T', '25', '-material', 'dna', '-multi', '-dangles', 'all', '-sodium', '0.05', '-magnesium', '0.0125', str(which)])

def computeMfe(strandList, name):
    global energyMap

    #hash to reduce reading/writing from files
    if name in energyMap:
        return energyMap[name]

    now = datetime.datetime.now()
    #create a file with input strands
    when = now.strftime('%b%d_%H%M%S')
    fname = DIR + name + '_' + when

    f = open(fname+'.in','w')
    f.write(str(len(strandList))+'\n');

    for i in range(len(strandList)):
        #print i,'#',strandList[i],'#\n'
        f.write(strandList[i]+'\n');

    for i in range(len(strandList)):
        f.write(str(i+1) + ' ')
    f.write('\n')
    f.close()

    #runmfe on it.
    runMfe(fname)

    #read the output file
    comment = '^%'
    emptyline = '^[ \t]*$'
    dotparen = '^[.()+]+$'
    impComment = '^% %'
    reComment = re.compile(comment)
    reEmpty = re.compile(emptyline)
    reDotParen = re.compile(dotparen)
    reImpComment = re.compile(impComment)
    ret = 0

    f = open(fname + '.mfe', 'r')
    lnum = 1

    #Checking each line in the file for the dot paren notation line
    while 1:
        line = f.readline()
        if line == "":
            break

        result = reImpComment.match(line)
        if result != None:
            line = f.readline()
            line = f.readline()
            ret = float(line)
            break

    #return the free energy value.
    f.close()
    energyMap[name] = ret
    return ret


def reverseComplement(strand):
    strand = strand.upper()
    strand = strand.replace('A','X')
    strand = strand.replace('T','A')
    strand = strand.replace('X','T')

    strand = strand.replace('G','X')
    strand = strand.replace('C','G')
    strand = strand.replace('X','C')
    strand = strand[::-1]
    return strand

#calculate energy and write into fenergy file
def calcAndStoreDomainEnergy(singleDomain, fenergy):

    if singleDomain in domainEnergyMap:
       return

    now = datetime.datetime.now()
    #create a file with input strands
    when = now.strftime('%b%d_%H%M%S')
    fname = DIR + singleDomain + '_' + when

    strand1 = PREFIX + singleDomain + SUFFIX
    strand2 = reverseComplement(strand1)
    both = [strand1,strand2]
    duplexMfe = computeMfe(both, strand1)

    energy = duplexMfe - computeMfe(BasicStrandList, BasicStrandName)

    f = open(fenergy,'a')
    f.write(singleDomain+' '+str(energy)+'\n');
    f.close()

    domainEnergyMap[singleDomain] = energy
#end calcAndStoreDomainEnergy

def computeAllDomainEnergies(current, toehold, fenergy):
    if(toehold < 1):
        calcAndStoreDomainEnergy(current, fenergy);
        return
    for j in range(4):
        curr_next = current+BASES[j]
        computeAllDomainEnergies(curr_next, toehold-1, fenergy)
    #end for
#end computeAllDomainEnergies

def main():
    if(len(sys.argv) != 3):
       print "Usage: python filename.py toeholdLength domainEnergy.txt"
       return
    readMappingFiles(sys.argv[2]);
    print "mapping files read"
    computeAllDomainEnergies('', int(sys.argv[1]), sys.argv[2])
    #endfor

#end main()

def readMappingFiles(fenergy):
    global domainEnergyMap
    emptyline = '^[ \t]*$'
    reEmpty = re.compile(emptyline)

    f = open(fenergy, 'r')
    for line in f:
        if(reEmpty.match(line)):
            continue
        w = line.split()
        domainEnergyMap[w[0]] = float(w[1])

    f.close()

#end readMappingFiles

#GlOBAL Declarations here
energyMap = dict()
domainEnergyMap = dict()
PREFIX='CATCGATCGATCG'
SUFFIX='CATCGATCGATCGATCG'
BasicStrandName=PREFIX+SUFFIX
BasicStrandList=[BasicStrandName,reverseComplement(BasicStrandName)]
DIR='runs/'
BASES='ATCG'

#Main function starts getting called from here
if __name__ == '__main__':
    main()
