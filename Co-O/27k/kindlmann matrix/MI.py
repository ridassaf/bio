import math

s = ''
for i in range(32):
    s+= '\t' + str(i)
s += '\n'

for i in range(32):
    s += str(i)
    for j in range(32):
        binary_i = bin(i)[2:]
        count_i = binary_i.count('1')
        prob_i = float(count_i)/5
        
        binary_j = bin(j)[2:]
        count_j = binary_j.count('1')
        prob_j = float(count_j)/5
        
        binary_i_j = bin(i&j)[2:]
        count_i_j = binary_i_j.count('1')
        prob_i_j = float(count_i_j)/5
        
        count_ni_j = count_j - count_i_j
        prob_ni_j = float(count_ni_j)/5
        
        count_i_nj = count_i - count_i_j
        prob_i_nj = float(count_i_nj)/5
        
        count_ni_nj = 5 - (count_i_j + count_ni_j + count_i_nj)
        prob_ni_nj = float(count_ni_nj)/5
        
        MI = 0
        if prob_i != 0 and prob_j != 0 and prob_i_j != 0:
            MI += prob_i_j*math.log(prob_i_j/(prob_i*prob_j),2)
        
        if prob_i != 1 and prob_j != 0 and prob_ni_j != 0:
            MI += prob_ni_j*math.log(prob_ni_j/((1.0-prob_i)*prob_j),2)
                
        if prob_i != 0 and prob_j != 1 and prob_i_nj != 0:
            MI += prob_i_nj*math.log(prob_i_nj/(prob_i*(1.0-prob_j)),2)
        
        if prob_i != 1 and prob_j != 1 and prob_ni_nj != 0:
            MI += prob_ni_nj*math.log(prob_ni_nj/((1.0-prob_i)*(1.0-prob_j)),2)
        
        s += '\t' + str(round(MI,2))
    s += '\n'

output = open('MI.txt', 'w')
output.write(s)
output.close()

