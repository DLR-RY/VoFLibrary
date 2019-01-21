from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="fileName", help="Open specified file")
args = parser.parse_args()
fileName = args.fileName

counterInitRes = 0
counterFinalRes = 0
counterIter = 0
avgInitRes = 0.0
avgFinalRes = 0.0
avgIter = 0.0

with open(fileName) as f:
    for line in f:
        if "final residual abs" in line:
            finalres = line.split()[4]
            avgFinalRes += float(finalres)
            counterFinalRes += 1
        	     #print (finalres)
        if "intial residual absolute" in line:
            initialres = line.split()[4]
            avgInitRes += float(initialres)
            counterInitRes += 1
        if "iterations" in line:
            initialres = line.split()[2]
            avgIter += int(initialres)
            counterIter += 1


print("avgFinalRes ",avgFinalRes/counterFinalRes)
print("avgIter ",avgIter/counterIter)
print("avgInitRes ",avgInitRes/counterInitRes)

