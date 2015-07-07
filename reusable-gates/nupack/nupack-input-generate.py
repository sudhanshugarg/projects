import re, string, sys
import subprocess

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
                subprocess.call(['mfe', '-T', '25', '-material', 'dna', '-multi', '-dangles', 'all', '-sodium', '0.05', '-magnesium', '0.0125', fname])


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

    def runMfe(self, which):
        subprocess.call(['mfe', '-T', '25', '-material', 'dna', '-multi', '-dangles', 'all', '-sodium', '0.05', '-magnesium', '0.0125', str(which)])

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
    ret = s.checkOutput(num)
    for i in range(num):
        print i, ret[i]

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

DIR='runs/'

inputStrands = [0 for i in range(100)]
for i in range(100):
    inputStrands[i] = dict()

#Main function starts getting called from here
if __name__ == '__main__':
    main()
