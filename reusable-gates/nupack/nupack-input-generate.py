import re, string, sys
import subprocess
import datetime

    
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

    def createInput(self,num):
        global inputStrands
        global nameMap
        global checkMap
        for whichInput in range(num):
            for i in range(len(checkMap)):
                w = checkMap[i].split()
                fname = DIR + str(whichInput) + '_' + w[0] + '_' + w[1]

                f = open(fname+'.in','w')
                f.write('2\n')
                f.write(inputStrands[whichInput][w[0]]+'\n');
                f.write(inputStrands[whichInput][w[1]]+'\n');
                f.write('1 2\n')
                f.close()
                runMfe(fname)
                #subprocess.call(['mfe', '-T', '25', '-material', 'dna', '-multi', '-dangles', 'all', '-sodium', '0.05', '-magnesium', '0.0125', fname])


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
                w = checkMap[i].split()
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

    def calculateDomainFreeEnergies(self):
        strandX = [0 for i in range(5)]
        strandY = [0 for i in range(5)]
        strandG = [0 for i in range(7)]
        #using complex X+G, and breaking up X into 2 strands.
        seqXG = [self.strand['X'], self.strand['G']]

        ldomain = len(self.strand['X'])/5
        Xl = self.strand['X'][:ldomain*3]
        Xr = self.strand['X'][ldomain*2:]

        seqXlG = [Xl, self.strand['G']]
        seqXrG = [Xr, self.strand['G']]

        #hash the info.
        strandX[1] = computeMfe(seqXG, 'XG') - computeMfe(seqXrG, 'XrG');
        strandX[3] = computeMfe(seqXG, 'XG') - computeMfe(seqXlG, 'XlG');

        #using complex Y+G, and breaking up Y into 2 strands.
        seqYG = [self.strand['Y'], self.strand['G']]

        Yl = self.strand['Y'][:ldomain*3]
        Yr = self.strand['Y'][ldomain*2:]

        seqYlG = [Yl, self.strand['G']]
        seqYrG = [Yr, self.strand['G']]

        strandY[1] = computeMfe(seqYG, 'YG') - computeMfe(seqYrG, 'YrG');
        strandY[3] = computeMfe(seqYG, 'YG') - computeMfe(seqYlG, 'YlG');

        #using complex Z+G, and breaking up Z into 2 strands.
        seqZG = [self.strand['Z'], self.strand['G']]
        Zl = self.strand['Z'][:ldomain*4]
        Zr = self.strand['Z'][ldomain*1:]

        seqZlG = [Zl, self.strand['G']]
        seqZrG = [Zr, self.strand['G']]

        strandX[2] = computeMfe(seqZG, 'ZG') - computeMfe(seqZrG, 'ZrG');
        strandY[2] = computeMfe(seqZG, 'ZG') - computeMfe(seqZlG, 'ZlG');

        seqXYG = [self.strand['X'], self.strand['Y'], self.strand['G']]
        seqDuplexG = [self.strand['G'], reverseComplement(self.strand['G'])]
        strandG[3] = computeMfe(seqXYG, 'GG') - computeMfe(seqDuplexG, 'XYG');

        ret = []
        ret.append(strandX)
        ret.append(strandY)
        ret.append(strandG)
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

def main():
    num = createInputFromNupackDesign('nd')
    #num = 1
    print "Number of Inputs: " + str(num)
    s = Strands()
    #print s.getStrand('MB')
    s.createInput(num)
    #ret = s.checkOutput(num)
    #for i in range(num):
    #    print i, ret[i]

            #inputStrands[input_num][nameMap[words[0]]] = words[1]
    s.strand['X'] = inputStrands[0]['X']
    s.strand['Y'] = inputStrands[0]['Y']
    s.strand['G'] = inputStrands[0]['G']
    s.strand['Z'] = inputStrands[0]['Z']
    energies = s.calculateDomainFreeEnergies()
    print energies
    #for i in range(3):
    #    for j in range(3):
    #        print energies[i][1], energies[i][2], energies[i][3], '\n'

#GlOBAL Declarations here
nameMap = dict()
nameMap['T4_X_strand'] = 'X'
nameMap['T4_Y_strand'] = 'Y'
nameMap['T4_Z_strand'] = 'Z'
nameMap['T4_G_strand'] = 'G'
nameMap['T4_GX_strand'] = 'GX'
nameMap['T4_GY_strand'] = 'GY'
nameMap['MB_strand'] = 'MB'

checkMap = dict()
checkMap[0] = 'MB X'
checkMap[1] = 'MB Y'
checkMap[2] = 'MB Z'
checkMap[3] = 'MB G'
checkMap[4] = 'MB GX'
checkMap[5] = 'MB GY'
checkMap[6] = 'Z X'
checkMap[7] = 'Z Y'
checkMap[8] = 'Z GX'
checkMap[9] = 'Z GY'

energyMap = dict()

DIR='runs/'

inputStrands = [0 for i in range(100)]
for i in range(100):
    inputStrands[i] = dict()

#Main function starts getting called from here
if __name__ == '__main__':
    main()
