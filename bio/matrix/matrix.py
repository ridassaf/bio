import sys
import re


input = open(sys.argv[1], 'r')

output = open('real-matrix.txt', 'w')

MI = {}
FIGS = []
for line in input:
    line = line.strip()
    words = line.split('\t')
    words[0] = words[0].replace('0','')
    words[1] = words[1].replace('0','')

    MI[(words[0], words[1])] = words[2]
    if words[0] not in FIGS:
        FIGS.append(words[0])
    if words[1] not in FIGS:
        FIGS.append(words[1])

for fig in FIGS:
    output.write('\t' + fig)
output.write('\n')
for fig in FIGS:
    output.write(fig)
    for fig2 in FIGS:
        if (fig,fig2) in MI:
            r = round(float(MI[(fig,fig2)]),2)
            output.write('\t' + str(r))
        else:
            output.write('\t-')
    output.write('\n')

output.close()
