#!/usr/local/bin/python

import collections
'''
import sys
inputFileName = sys.argv[1]
outputFileName = sys.argv[2]
'''

################## Import ##################

# Open input file
# inputFileName = input('Enter input file name:')
try:
    # f= open(inputFileName, 'r')
    # f = open('dna.example.fasta', 'r')
    f = open('dna1.fasta', 'r')
except IOError:
    print('The file does not exist!')
# f.seek(0)

# Read sequences into a dictionary
seqs = {}
for line in f:
    line = line.rstrip()
    if line[0] == '>':
        header = line.split()
        name = header[0][1:]
        seqs[name] = ''
    else:
        seqs[name] = seqs[name] + line

# Print the total # of sequences.
print('There are', len(seqs), 'sequences.')


################## Longest sequence ##################

# Find the longest length and print the identifier(s).
print('\n')
maxLen, maxKey = max((len(v), k) for k, v in seqs.items())
print('The longest sequence is:', maxLen, 'nt')
for k, v in seqs.items():
    if len(v) == maxLen:
        print(k, len(v))

# Find the shortest length and print the identifier(s).
print('\n')
minLen, minKey = min((len(v), k) for k, v in seqs.items())
print('The shortest sequence is:', minLen, 'nt')
for k, v in seqs.items():
    if len(v) == minLen:
        print(k, len(v))


################## ORF ##################

# Define start and stop codons.
print('\n')
startCodon = {'ATG'}
stopCodon = {'TAA', 'TAG', 'TGA'}

#p = 1
seq = 'ATGATGCCCTAGcca'
seqlist = ['ATGATGCCCTAGcca', 'aATGacagctTAGATGCC', 'aATGaATGCCCTAG']
seqdic = collections.OrderedDict({'S1':'ATGATGCCCTAGcca', 'S2':'aATGacagctTAGATGCC', 'S3':'aATGaATGCCCTAG'})
testSeq = 'ATGCGCGCGTCCGCTCGCCGCGATAGATCGACTGCACGATGCCGGTCGTGTCGAGCCCTTCGAGCAGCGCGCGGCGCACCTTCACGTCCTGCGTGGGGCCGTAGCCGACGTTCAGGAACAGCGAGTACGGCGTGCCGGTGTTCAGCGCGTGCTGGTAGGTGAAGTCCGGGTTCTTGCGGATCAATGCCGCGTCGTTGCCGGATACGCCCTCGATCACGTCGACCTGGCCGGACAGCAGCGCGCCGGTGCGCACCGACGATTCCGGCAGGAACCGGTAGACGATCCGGTCGAGGTAGGCCGGGCCTTGATGGCCCGCCGTCGGCGGCGCCCAGCGATAGTCGGGATTCTTCGCGAATTCGATTTCCTGCCCCTTCACGTAGCGGCGCAGGATGAACGGGCCGGTGCCGGCGACGTCGGTGCCGCCCGCTTTCACGCTCGGCCGGTCGAACGCGTTCGGCGACAGCAGCGCGAACCGGCCGGCGAAGGTCAGGAATGGCGCATACGGCGCGTCGAGCGTGATCACGACCGTGCGCGCGTCCGGCGTCTTCACGTCCGCGACGTGCAGGATCTGCGCCGCGAGCTGGCTGCCGCCCGAATACGCCGGATTGCGTGCGTTCAGGAAGTTGCGCGCGACGGCCTGCGCGTCGAGCGGCGCGCCGTCGCTGAACTTCACGCCGTCGCGCAGCGTGAACGTGTAGGTCTTGCCGTCGGGCGACACCTTGTAGCCGGTCGCGAGCCACGGCACGAAGCTGCCGTCCGGCTTGCGCGCGAGCAGGCTCTCATACACGTTGCGCAGCAGCAGCTCGACCTTGTCCTGGCCGTTCAATTGCGGATTCAGCGTGCGCGGCTCGGATTCGACGCCCCAGGTCAGCGTGCCGCCCGTGACGGGCGTGCCCTGGGCGAACGCCGCGGGCGGCAAGGCGAGGGCGAACAGCAGCGCGGCCGATGCCGTGCGGCGCGACGCGCGGAGCGTGAAATCCGGAAAAGCGAAGGTCATGCGGGGCAATCCTTTACAAGTGGCATGGCGGGTGGAACGGGGCACGACGAAGCGCGGTGCGAGAGGCGTGTCGAAACGGAGCGCGTGGAACGTCATGCGCTCGCGCGAACACGCCCGGCGCCCTGGGCGGCATCCGCCGCCGGCGGTTCCCGGTAGAGATCCGGGATCAGGATCTTTCTCAGCCAGACGAGCGGCGTGCCGAGCGTGTTCGTCGAATGGTCGAACGCGACATAGCCTTGCCGCTCGTAGTACGGCTGCAGCCACGGATGATTGGTGGCGGTCGCGAGATAGGTCGCGGGCGCCTTGACCTGTTCGGTCAGGAATTGTGCTTCGACGTGGTCGAGCAGTGCGAGGCCGAGTCCCTGCCGCTTGAACGCGGGCGCCACCGCGAACCAGTGCAGGAACGGATACGGCGTGCGGTGGCGCTCGCCGGACACCCACGGGAACCGCACGGTCAGCGTCGCGGCAGGCGCGCGGCGGCCGTCGGCGCCGGCGCGTTCGAGCACGAAGGTCGTTTCGCGCGCAACGGTCTGCGCGACGTGCTCGATCGGCGAGCGAACGATGTTGAAGAGCACGCCCTGCGGCACGAGCGGCGCATACGCGTCCTGCAGCAACGCATAAAGGGCCGGCACGTCGTCGGCGCCGGCAATGCGGATCGTATCGGTCATCGTGCGGCGCCTTCGAGGTCGGGTTGACGGGTGACGGCGGCGGTGTGCGGTGTCCGTGGCGTTTCTTCGGTCGCGCGTCGTGACGCCGCGGGGTCGCCGCTGTCCGCCGACGCACCGGCCAGAAGCCGCAGCGACGCAATCCTCGGCGCTGCATCGGCAAGCGGTGTGTCGACGACGAATTCCGCGATCCCGTGGCGGCGATGCAGCGTGTCGAGTTCGCGCAGCACGTCGTCGCGCGTGCCGTGCACGATCCGCGTCTCGCGCGCGTCGATCCGATGCGGCCCGTCCCCCGCCTGGCGGATGTACGCGTGTGCCTGTTCGAGACTGCCGACGTTGACCGCCGGCGTGTCGCCCACTTGCACGCGAAAGCCCTGATGCCCGCTCGCGAGTCGCGCGGCCGCTTCGCGCGTGGCTGCGACGATCACGGATACGGCGAGCAGCGGCGTGCGGCCGCCCGATGCGTGGCGGTATGCGTCGAAGCTCCGTTCCAGATCGGCGGTGGATGCGTTGAGGTGCGACGCATAGACGAACGCCCAGCCGAGCTGGGCGGCCAGCGCGGCGCTTTCGGCGCTCGCGCCGAGCAGGAACGGTTGGGCGGCCACGGGCGGAAACGGTGTCGCGTGCAGCGTGGCGTCCGTGTGCCCCGCGTCGTCCGGCGGCGTGCCGCTGCCGTCGAGATGGCGCTTCAGCTCGGCTGCGAGTTCGGCGAACGACGCGCGTTGCGCCGGGTCGTAGCCGGCCTGCAGCGCGCGCGTCGACAACGGCAGCCCGCCCGGCGCCTTGCCGATGCCGAGATCGACGCGATCGGGCGCAAGTGACGCGAGCAGGTTGAAGTTCTCCGCGATCTTGTACGGGCTGTAGTGCTGCAGCATCACGCCGCCCGATCCGACACGAATTCGCCGCGTGTGTGCGAGCACATGCGCGATCAGGATTTCCGGCGACGGGCATGCGAGCGCGGCCGTGTTGTGATGCTCGGCGAGCCAGTAGCGGTGATAGCCCGCATCGTCGGCGAAGCGCGCGAGCGACACGCTGCGCGCGAGTGCCTCCACGGCGGTTTCGCCCGGTGCGACGGGGCTCTTGTCGAGCACGCTCAGGCGAAACGGGAGGGCGGGTAGCGCTGCGGCGGAATGAGGGGCGTCGGACGTCATGAGGGTGCGACCGGCGCGGTTCGAAGCGGCACGATCAACTGGTTAGGCATGCCGGATTATTTTTCGACGCCGTTTGCACGGTCAAACAATCAAGAAGAATATCGAAAGCCGAAAATATTCGTTCCGGCTCACCGGGTGGCCTGCGTTATGGTCGATGGCTACGCGCCTTGCGGCGCACGCATGGACCCTGCCGGATGCCCCGGCGCGTGATCGCCGGCAGCATGGTGCCGGCTCGGGTCCGCCGATCGAACGCAATCGTGGGGAATACCGATGAACGAGCGAATCATGTCAACCGGCCGACCTGCGGCGCATGCCGGAGCCGCATCATGAGCGATACACGCAGACGGCTGAAGACGCTCGTCGGCGCCGCGAGCGTGCTGGCTGCGCCGGGCGATGCGGATCCTGCGCAAGGCATGCGTCGCGCGGCGGCGCTCGCGCGCAAGGCGGAAGCCGAAAAGATCACCGGCCTCTTCACGGCCGACCTGTTGCAGGCCGATCCGGCCGGGCTCGCCGGCAGCACCGGCAGCCAGGAGCCGATCGTCGCGCTGGCCGCGCTGAGCCAGGCGACGTCGGCCATCGGGCTCGTCGCGACGGTGTCGACCACGTATCACCATCCGTACAACCTTGCGCGGCTGATCGGCACGCTCGATCACGTGAGCGGCGGCCGCGCCGCGTGGAACGCGGTGACGTCGTCGGTCGGCGAGGAGAATTTCGGCGACGCCGCGTTGCCCGACCCGGAGCAGCGCTATGCGCGCGCCGCCGAATTCGTCGAAGTCGTCAACGCGCTGTTCGACGCGAACGATCCCGACGCGGTGCGACGCACGCCGAGCGGCGGCGTGTCGGTCGATCCGCTCAAGCTCGGGCCGATCGATTACCGTGGCGCGCACTTCCAGGTGCAGGGGCCGTTGAACGTGCCGCCTTCGCCGCAACGGCGCCCCGTGCTGTTTCAGGCGGGGCAGTCGGCGAGCGGCGTGACGCTCGGCGCGCGGTTTGCGGAAGTCGTCTACACGTCGCAGCCGACGCTCGAGGATGCGCGCACGTTCGTCGGCGAACTGCACCGGCAGGCGGCCGCATTCGGCCGGGCAGGGCGGCTGCCGCTCGTGATGAACTCGTTCCACTCGGTGATCGGCGAATCCGACGCCGACGTCGCGCGACGGCTGCGCGACAAGCATGA'

# :::::::: Define Function: Find ORFs in a given RF in one sequence ::::::::

def findORF(seq, p):
    ORFnum = 0
    ORF = collections.OrderedDict()
    for nt in range(p - 1, len(seq), 3):
        if seq[nt:nt + 3] in startCodon:
            #print('Find ATG at position', nt)
            noStop = 1
            for nt2 in range(nt + 3, len(seq), 3):
                if (seq[nt2:nt2 + 3] in stopCodon) * noStop:
                    noStop = 0
                    #print('find a ORF:', seq[nt:nt2+3])
                    ORF[str(nt+1)] = seq[nt:nt2 + 3]
                    ORFnum += 1
                # print('Find', ORFnum, 'ORF(s) in reading frame', p)
                # print(ORF)
    return (ORF)

#print('******************')
#testSeqORF = collections.OrderedDict()
#testSeqORF = findORF(testSeq, 1)
#print('******************')


# :::::::: Define Function: Find ORFs in all RFs in one sequence ::::::::

def findORFinAllRF(seq, rf_flag = 'all'):
    ORFinAllRF = collections.OrderedDict()
    if rf_flag == 'all':
        for p in range(1, 4):
            rf = 'RF'+str(p)
            ORFinAllRF[rf] = findORF(seq, p)
        return (ORFinAllRF)
    else:
        rf = 'RF'+ str(rf_flag)
        ORFinAllRF[rf] = findORF(seq, int(rf_flag))
        return (ORFinAllRF)

ORFinAllRF = findORFinAllRF(seq)
# print('\nORFs in all reading frame:',ORFinAllRF)

# :::::::: Define Function: Find ORFs in all RFs in all sequences ::::::::

def findORFinAllRFinFile(seqdic, rf_flag):
    allORF = collections.OrderedDict()
    for k in seqdic.keys():
        allORF[k] = findORFinAllRF(seqdic[k], rf_flag)
    return (allORF)


# :::::::: Call Function: Find ORFs in all RFs in one sequence ::::::::

#allORF = collections.OrderedDict()
#allORF = findORFinAllRFinFile(seqs, 'all')

# :::::::: Print OFR result ::::::::

'''
for k in allORF.keys():
    print(k)
    for r in allORF[k].keys():
        print(r, ':', allORF[k][r])
    print('\n')
'''

def flatten(d, parent_key='', sep='_'):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def maxORF(allORF):
    allORFflat = flatten(allORF)
    maxORFLen, maxORFSeq, maxORFKey = max((len(v), v, k) for k, v in allORFflat.items())
    #print(maxORFKey, maxORFLen)
    print('The longest ORF has', maxORFLen, 'nt.')
    for k, v in allORFflat.items():
        if len(v) == maxORFLen:
            print('Identifier:', k, '\tLength:', len(v))
            print(v)
    print('\n')


allORF = collections.OrderedDict()
allORF = findORFinAllRFinFile(seqs, 'all')

print('\n######### Longest ORF in all sequences. #########')
maxORF(allORF)

print('\n######### Longest ORF in each sequence. #########')
for k in allORF.keys():
    print('For sequence', k)
    allORFflat = flatten(allORF[k])
    if len(allORFflat) != 0:
        maxORFLen, maxORFSeq, maxORFKey = max((len(v), v, k) for k, v in allORFflat.items())
        #print(maxORFKey, maxORFLen)
        print('The longest ORF has', maxORFLen, 'nt.')
        for k, v in allORFflat.items():
            if len(v) == maxORFLen:
                print('Identifier:', k, '\tLength:', len(v))
    else:
        print('No ORF is founded.')
    print('\n')


print('\n######### Longest ORF in each reading frame. #########')
# Extract all ORF of a given reading frame in all seqeunces.
def ORFFrom1RF(seqdic, rf_flag):
    seqsFrom1RFdic = collections.OrderedDict()
    rf = 'RF' + rf_flag
    for s in seqdic.keys():
        seqsFrom1RFdic[s] = seqdic[s][rf]
    return(seqsFrom1RFdic)

def maxORFFrom1RF(allORF, rf_flag):
    print('For reading frame', rf_flag, 'in all sequences.')
    seqsFrom1RF = collections.OrderedDict()
    seqsFrom1RF = ORFFrom1RF(allORF, rf_flag)
    maxORF(seqsFrom1RF)

maxORFFrom1RF(allORF, '3')


################## Motif ##################

def countMotif(seq, motif):
    count = start = 0
    while True:
        start = seq.find(motif, start) + 1
        if start > 0:
            count+=1
        else:
            return count

def findMotif(seqdic, motifLen):
    motif = collections.OrderedDict()
    for s in seqdic.keys():
        for nt in range(0, (len(seqdic[s]) - motifLen)):
            m = seqdic[s][nt:(nt + motifLen)]
            if m not in motif.keys():
                num = 0
                for k in seqdic.keys():
                    seqStr = str(seqdic[k])
                    #num += seqStr.count(m)
                    num += countMotif(seqStr, m)
                if num > 1: motif[m] = num
    return(motif)

motifTest = {'1':'acacacacadadadad', '2':'adadadadacacacac', '3':'efgefgefg'}
motifRes = findMotif(seqs, 12)
print(motifRes)



maxMotifValue, maxMotifKey = max((v, k) for k, v in motifRes.items())
print(maxMotifKey, maxMotifValue)
for k, v in motifRes.items():
    if v == maxMotifValue:
        print('Motif:', k, '\tOccurence:', v)

################# Output ##################

# Close input file
f.close()
