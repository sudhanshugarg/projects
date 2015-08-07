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
    when = now.strftime('%b%d_%H%M')
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


class Strands(object):
    def __init__(self):
        self.strand = {
            'X': '',
            'Y': '',
            'Z': '',
            'G': '',
            'GX': '',
            'GY': '',
            'MB': ''
        }

    def printme(self):
        greeting = ', '.join("%s=%r" % (key,val) for (key,val) in k.iteritems())
        print greeting

    def getStrand(self, which):
        if which != 'MB':
            return self.strand[which]
        else:
            return 'MB:'+self.strand[which]

    def sieveBadInput(self,num):
        global inputStrands
        global nameMap
        global checkMap
        actualNum = 0
        for whichInput in range(num):
            #check if all domains match required criteria.
            self.strand['X'] = inputStrands[whichInput]['X']
            self.strand['Y'] = inputStrands[whichInput]['Y']
            self.strand['G'] = inputStrands[whichInput]['G']
            self.strand['Z'] = inputStrands[whichInput]['Z']
            energies = self.calculateDomainFreeEnergies(str(whichInput))
            print "Input: ", whichInput, " energies: ", energies
            if energies[3] != 7:
                print "This input, ", whichInput, ", does not work: ", energies[3], "\n"
                #continue


            print "This input, ", whichInput, ", Works!\n"

            for i in range(len(checkMap)):
                w = checkMap[i].split(',')
                fname = DIR + str(actualNum) + '_' + w[0] + '_' + w[1]

                f = open(fname+'.in','w')
                f.write('2\n')
                f.write(inputStrands[whichInput][w[0]]+'\n');
                f.write(inputStrands[whichInput][w[1]]+'\n');
                f.write('1 2\n')
                f.close()
                runMfe(fname)
                #subprocess.call(['mfe', '-T', '25', '-material', 'dna', '-multi', '-dangles', 'all', '-sodium', '0.05', '-magnesium', '0.0125', fname])

            actualNum += 1

        return actualNum

    def checkOutput(self,num):
        comment = '^%'
        emptyline = '^[ \t]*$'
        dotparen = '^[.()+]+$'
        reComment = re.compile(comment)
        reEmpty = re.compile(emptyline)
        reDotParen = re.compile(dotparen)
        ret = []

        for whichInput in range(num):
            ret.append(0)
            for i in range(len(checkMap)):
                w = checkMap[i].split(',')
                fname = DIR + str(whichInput) + '_' + w[0] + '_' + w[1] + '.mfe'

                flag = 0
                f = open(fname, 'r')
                lnum = 1

                #Checking each line in the file for the dot paren notation line
                for line in f:
                    result = reDotParen.match(line)
                    if result != None:
                        numparen = line.count('(')
                        if numparen > 5:
                            flag = 1
                        #As soon as we find the dot paren line, we process it and
                        #then break out of the loop
                        break
                        #sys.stdout.write (str(lnum) + ':' + line)

                #Done reading the file
                if flag == 1:
                    print "Input: " + str(whichInput) + " : " + fname + " DOESNT WORK"
                else:
                    print "Input: " + str(whichInput) + " : " + fname + " works"
                    ret[-1] += 1

                #Now, close the file and continue with the next file in the same input
                f.close()

        return ret

    def calculateDomainFreeEnergies(self, inputName):
        strandX = [0 for i in range(5)]
        strandY = [0 for i in range(5)]
        strandG = [0 for i in range(7)]
        #using complex X+G, and breaking up X into 2 strands.
        seqXG = [self.strand['X'], self.strand['G']]

        ldomain = len(self.strand['X'])/5
        Xl = self.strand['X'][:ldomain*3]
        Xl4 = self.strand['X'][:ldomain*4]
        Xr = self.strand['X'][ldomain*2:]
        Xr4 = self.strand['X'][ldomain:]

        seqXlG = [Xl, self.strand['G']]
        seqXrG = [Xr, self.strand['G']]
        seqXl4Xl4 = [Xl4, reverseComplement(Xl4)]
        seqXr4Xr4 = [Xr4, reverseComplement(Xr4)]
        seqXX = [self.strand['X'], reverseComplement(self.strand['X'])]

        #hash the info.
        strandX[0] = round(computeMfe(seqXX, 'XX'+inputName) - computeMfe(seqXr4Xr4, 'Xr4'+inputName),2);
        strandX[1] = round(computeMfe(seqXG, 'XG'+inputName) - computeMfe(seqXrG, 'XrG'+inputName),2);
        strandX[3] = round(computeMfe(seqXG, 'XG'+inputName) - computeMfe(seqXlG, 'XlG'+inputName),2);
        strandX[4] = round(computeMfe(seqXX, 'XX'+inputName) - computeMfe(seqXl4Xl4, 'Xl4'+inputName),2);

        #using complex Y+G, and breaking up Y into 2 strands.
        seqYG = [self.strand['Y'], self.strand['G']]

        Yl = self.strand['Y'][:ldomain*3]
        Yl4 = self.strand['Y'][:ldomain*4]
        Yr = self.strand['Y'][ldomain*2:]
        Yr4 = self.strand['Y'][ldomain:]

        seqYlG = [Yl, self.strand['G']]
        seqYrG = [Yr, self.strand['G']]
        seqYl4Yl4 = [Yl4, reverseComplement(Yl4)]
        seqYr4Yr4 = [Yr4, reverseComplement(Yr4)]
        seqYY = [self.strand['Y'], reverseComplement(self.strand['Y'])]

        #hash the info.
        strandY[0] = round(computeMfe(seqYY, 'YY'+inputName) - computeMfe(seqYr4Yr4, 'Yr4'+inputName),2);
        strandY[1] = round(computeMfe(seqYG, 'YG'+inputName) - computeMfe(seqYrG, 'YrG'+inputName),2);
        strandY[3] = round(computeMfe(seqYG, 'YG'+inputName) - computeMfe(seqYlG, 'YlG'+inputName),2);
        strandY[4] = round(computeMfe(seqYY, 'YY'+inputName) - computeMfe(seqYl4Yl4, 'Yl4'+inputName),2);

        #using complex Z+G, and breaking up Z into 2 strands.
        seqZG = [self.strand['Z'], self.strand['G']]
        Zl = self.strand['Z'][:ldomain*4]
        Zr = self.strand['Z'][ldomain*1:]

        seqZlG = [Zl, self.strand['G']]
        seqZrG = [Zr, self.strand['G']]

        strandX[2] = round(computeMfe(seqZG, 'ZG'+inputName) - computeMfe(seqZrG, 'ZrG'+inputName),2);
        strandY[2] = round(computeMfe(seqZG, 'ZG'+inputName) - computeMfe(seqZlG, 'ZlG'+inputName),2);

        seqXYG = [self.strand['X'], self.strand['Y'], self.strand['G']]
        seqDuplexG = [self.strand['G'], reverseComplement(self.strand['G'])]
        strandG[3] = round(computeMfe(seqDuplexG, 'GG'+inputName) - computeMfe(seqXYG, 'XYG'+inputName),2);

        #lower and higher thresholds
        lowT = -7
        highT = -9
        goodDomains = 0
        for i in range(3):
            j = i+1
            if(strandX[j] <= lowT and strandX[j] >= highT):
                goodDomains += 1
            if(strandY[j] <= lowT and strandY[j] >= highT):
                goodDomains += 1

        if(strandG[3] <= lowT and strandG[3] >= highT):
            goodDomains += 1

        ret = []
        ret.append(strandX)
        ret.append(strandY)
        ret.append(strandG)
        ret.append(goodDomains)
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

def createInputFromNupackDesign(fname):
    global inputStrands
    global nameMap
    f = open(fname, 'r')
    emptyline = '^[ \t]*$'
    defect = 'ensemble defect'
    reEmpty = re.compile(emptyline)
    reDefect = re.compile(defect)

    input_num = -1
    for line in f:
        result = reDefect.search(line)
        if result != None:
            input_num += 1
            continue
        result = reEmpty.match(line)
        if result == None:
            words = line.split()
            #print "hi:" + words[0] + ":" + words[1]
            inputStrands[input_num][nameMap[words[0]]] = words[1]
    f.close()
    return input_num+1
#end createInputFromNupackDesign()

def createInputFromCurlToNupack(fname):
    '''
    #read nupack script from this file.
    with open (fname, 'r') as content_file:
        design = content_file.read()
    '''

    multipart = {
                'design_job[target_structure]': open(fname,'rb')
            }
    payload = {
            'preview_token' : '',
            'design_job[nucleic_acid_type]' : 'DNA',
            'design_job[temperature]' : '25.0',
            'design_job[number_of_trials]' : '10',
            'design_job[rna_parameter_file]' : 'rna1995',
            'design_job[dna_parameter_file]' : 'dna1998',
            'design_job[dangle_level]' : '2',
            'design_job[na_salt]' : '0.5',
            'design_job[mg_salt]' : '0.0125',
            'design_job[dotplot_target]' : '0',
            'design_job[prevented_strings]' : 'SSSS\nAAAAA\nAAAA\nGGG',
            'design_job[email_address]' : 'sgarg@cs.duke.edu',
            'commit' : 'Design'
    }
    hdrs = {
            '''
            'Cache-Control':'max-age=0',
            'Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
            'Origin':'null',
            'User-Agent':'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36',
            'Content-Type':'multipart/form-data; boundary=----WebKitFormBoundaryRvFIpcbpreh4yDoF',
            'DNT':'1',
            'Accept-Encoding':'gzip, deflate',
            'Accept-Language':'en-US,en;q=0.8',
            '''
    }
    submit_page1 = 'http://www.nupack.org/design/new'
    rd = requests.post(submit_page1, data=payload, headers=hdrs, files=multipart)
    print rd.text.encode('utf-8')

    #extract the return url from this return data.

    #keep trying the return url until the design is complete.
    flag = 1
    #sleep for 10 seconds.
#end createInputFromCurlToNupack()

def sieveEachStrand(num, fenergy):
    global inputStrands
    global nameMap
    global checkMap
    global domainEnergyMap

    reOutput = '_Z$'
    reGate = '_G$'
    reNotGate = '_NG$'
    reStrandisOutput = re.compile(reOutput)
    reStrandisGate = re.compile(reGate)
    reStrandisNOTGate = re.compile(reNotGate)

    for i in range(num):
       #check if all domains match required criteria.
       for key in inputStrands[i]:
          print "current key:", key
          if(reStrandisOutput.search(key)):
             print "matched output"
             storeDomainFreeEnergy(inputStrands[i][key], 5, fenergy)
          if(reStrandisGate.search(key)):
             print "matched gate"
             storeDomainFreeEnergy(inputStrands[i][key], 7, fenergy)
          if(reStrandisNOTGate.search(key)):
             print "matched not gate"
             storeDomainFreeEnergy(inputStrands[i][key], 4, fenergy)

       #endfor
    #endfor
#end sieveEachStrand

def storeDomainFreeEnergy(strand, domains, fenergy):
   global domainEnergyMap
   domainLength = len(strand)/domains
   for i in range(domains):
      singleDomain = strand[i*domainLength:(i+1)*domainLength]
      calcAndStoreDomainEnergy(singleDomain, fenergy)
   #endfor
#end storeDomainFreeEnergy

#bound the domain by two pre-defined domains, and
#calculate the difference in free energies, and write
#this into the file
#calculate energy and write into fenergy file
def calcAndStoreDomainEnergy(singleDomain, fenergy):

    if singleDomain in domainEnergyMap:
       return

    now = datetime.datetime.now()
    #create a file with input strands
    when = now.strftime('%b%d_%H%M')
    fname = DIR + singleDomain + '_' + when

    strand1 = PREFIX + singleDomain + SUFFIX
    strand2 = reverseComplement(strand1)
    both = [strand1,strand2]
    duplexMfe = computeMfe(both, strand1)

    energy = duplexMfe - computeMfe(BasicStrandList, BasicStrandName)

    f = open(fenergy,'a')
    f.write(singleDomain+':'+str(energy)+'\n');
    f.close()

    domainEnergyMap[singleDomain] = energy
#end calcAndStoreDomainEnergy

def createPairedInput(num):
   global inputStrands
   global nameMap
   global checkMap
   actualNum = 0
   for whichInput in range(num):
      for i in range(len(checkMap)):
         w = checkMap[i].split(',')
         fname = DIR + str(actualNum) + '_' + w[0] + '_' + w[1]

         f = open(fname+'.in','w')
         f.write('2\n')
         f.write(inputStrands[whichInput][w[0]]+'\n');
         f.write(inputStrands[whichInput][w[1]]+'\n');
         f.write('1 2\n')
         f.close()
         runMfe(fname)
   actualNum += 1
   return actualNum

def main():
    if(len(sys.argv) != 4):
       print "Usage: python filename.py sequence.np xor.map domainEnergy.txt"
       return
    readMappingFiles(sys.argv[2], sys.argv[3]);
    print "mapping files read"
    num = createInputFromNupackDesign(sys.argv[1])
    print "Number of Inputs: " + str(num)
    s = Strands()
    #print s.getStrand('MB')
    sieveEachStrand(num, sys.argv[3])
    print "sieving strands done"
    createPairedInput(num)
    ret = s.checkOutput(num)
    for i in range(num):
        print i, ret[i]

#end main()

def readMappingFiles(fmap, fenergy):
    global nameMap
    global checkMap
    global domainEnergyMap
    f = open(fmap,'r')

    emptyline = '^[ \t]*$'
    reEmpty = re.compile(emptyline)

    checkPattern = '%checkMapPattern_Here%'
    reCheckPattern = re.compile(checkPattern)
    flag = 0
    count = 0

    for line in f:
        if(reEmpty.match(line)):
            continue

        if(reCheckPattern.search(line)):
            flag = 1
            continue

        w = line.split()
        if (flag == 0):
            nameMap[w[0]] = w[1]
        else:
            checkMap[count] = w[0]
            count += 1

    f.close()

    f = open(fenergy, 'r')
    for line in f:
        if(reEmpty.match(line)):
            continue
        w = line.split(':')
        domainEnergyMap[w[0]] = float(w[1])

    f.close()

#end readMappingFiles

#GlOBAL Declarations here
nameMap = dict()
checkMap = dict()
energyMap = dict()
domainEnergyMap = dict()
PREFIX='CATCGATCGATCG'
SUFFIX='CATCGATCGATCGATCG'
BasicStrandName=PREFIX+SUFFIX
BasicStrandList=[BasicStrandName,reverseComplement(BasicStrandName)]
DIR='runs/'

inputStrands = [0 for i in range(100)]
for i in range(100):
    inputStrands[i] = dict()

#Main function starts getting called from here
if __name__ == '__main__':
    main()
