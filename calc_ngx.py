from optparse import OptionParser
import sys
from joblib import Parallel, delayed
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-p","--prefix",action = 'store',type = 'string',dest = 'prefix',help = "")
parser.add_option("-t","--threadN",action = 'store',type = 'int',dest = 'threadN',help = "")
parser.add_option("-g","--genomeSize",action = 'store',type = 'int',dest = 'genomeSize',help = "")
(opt, args) = parser.parse_args()
if opt.threadN == None or opt.prefix == None:
    print("Error! usage: calc_ngx.py -t 24 -p keumgang -g 16000000000(option)")
    sys.exit()
batchN = opt.threadN
prefix = opt.prefix
genomeSize = opt.genomeSize
#####################################################################################
class MFA:
    def __init__(self, fileName):
        self.seqName_LIST = []
        self.seq_DICT = {}

        fin = open(fileName)
        for line in fin:
            if line.startswith('>') == True:
                seqName = line.rstrip('\n').split('\t')[0].split(' ')[0][1:]
                self.seqName_LIST += [seqName]
                self.seq_DICT[seqName] = []
            else:
                sequence = line.rstrip('\n').upper()
                self.seq_DICT[seqName] += [sequence]
        fin.close()

        for seqName in self.seqName_LIST:
            self.seq_DICT[seqName] = ''.join(self.seq_DICT[seqName])
        
        print('[done] mfa file read - ' + fileName, flush=True)

#####################################################################################
def run_batch(seq_DICT, batchN):

    def run_single(seqName, sequence):
        result_LIST = []
        flag = False
        for pos, nucl in enumerate(sequence):
            if nucl != 'N':
                if flag == False:
                    flag = True
                    sPos = pos
                    ePos = pos
                else:
                    ePos = pos
            else:
                if flag == False:
                    pass
                else:
                    flag = False
                    result_LIST += [(seqName, ePos - sPos + 1, sPos, ePos)]
        if flag == True:
            result_LIST += [(seqName, ePos - sPos + 1, sPos, ePos)]
            
        print('[done] find contig - ' + seqName, flush=True)
        return (seqName, result_LIST)
    
    batchResult = Parallel(n_jobs=batchN)(delayed(run_single)(seqName, sequence) for seqName, sequence in seq_DICT.items())

    result_DICT = {}
    for seqName, result_LIST in batchResult:
        result_DICT[seqName] = result_LIST

    return result_DICT

#####################################################################################
class Statistics_Length:
    def __init__(self, length_LIST):
        self.length_LIST = length_LIST

        self.sortedLength_LIST = sorted(self.length_LIST, reverse=True)
        self.totalLength = sum(self.sortedLength_LIST)

        self.maxLength = self.sortedLength_LIST[0]
        self.minLength = self.sortedLength_LIST[-1]

    def get_total(self):
        return self.totalLength
    
    def get_lengthN(self):
        return len(self.sortedLength_LIST)

    def get_max(self):
        return self.maxLength
    
    def get_min(self):
        return self.minLength

    def get_N(self, x):
        sumLength = 0
        for idx, length in enumerate(self.sortedLength_LIST):
            sumLength += length
            if float(sumLength) / self.totalLength >= (float(x)/100):
                return length, idx+1
    
    def get_NG(self, x, genomeSize):
        sumLength = 0
        for idx, length in enumerate(self.sortedLength_LIST):
            sumLength += length
            if float(sumLength) / genomeSize >= (float(x)/100):
                return length, idx+1
        return length, idx+1

#####################################################################################
ref = MFA(prefix + '.fa')
fout = open(prefix + '.length', 'w')
length = 0
for seqName, sequence in ref.seq_DICT.items():
    length += len(sequence)
fout.write(str(length) + '\n')
fout.close()

result_DICT = run_batch(ref.seq_DICT, batchN)

length_LIST = []
for seqName in ref.seqName_LIST:
    for result in result_DICT[seqName]:
        length_LIST += [result[1]]

#Statistics_Length
statistics_Length = Statistics_Length(length_LIST)

fout = open(prefix + '.contig', 'w')
for seqName in ref.seqName_LIST:
    for result in result_DICT[seqName]:
        fout.write('\t'.join(map(str, result)) + '\n')
fout.close()


#NG1 ~ NG100
fout = open(prefix + '.' + 'N', 'w')
for x in range(1, 101):
    N, L = statistics_Length.get_N(x)
    fout.write(f'N{x}\t{N}\tL{x}\t{L}\n')
fout.close()


if opt.genomeSize != None:
    #NG1 ~ NG100
    fout = open(prefix + '.' + 'NG', 'w')
    for x in range(1, 101):
        NG, LG = statistics_Length.get_NG(x, genomeSize)
        fout.write(f'NG{x}\t{NG}\tLG{x}\t{LG}\n')
    fout.close()



