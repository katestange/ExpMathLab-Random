import math
from scipy import stats
import itertools

significance = 0.05
significance2 = 0.01

def stringify(l):
    a = ""
    for c in l:
        a += str(c)
    return a

def to_int(b):
    total = 0
    i = 0
    for bit in reversed(b):
        total += bit*(2**i)
        i+=1
    return total

def tobin(seq):
    a = []
    for b in seq:
        t = format(b, "b")
        for c in t:
            a.append(int(c))
    print(len(a))
    return a



def serial(seq):
    if len(seq) < 21:
        return "error, sequence too small."
    n1 = seq.count(1)
    n0 = seq.count(0)
    n00 = 0
    n01 = 0
    n10 = 0
    n11 = 0
    for i in range(len(seq)-1):
        subseq = str(seq[i])+str(seq[i+1])
        if (subseq == "00"):
            n00 += 1
        elif(subseq == "01"):
            n01 += 1
        elif(subseq == "10"):
            n10 += 1
        elif(subseq == "11"):
            n11 += 1
    statistic = ((4/(len(seq)-1.0))*((n00**2)+(n01**2)+(n10**2)+(n11**2)))-((2*((n0**2)+(n1**2)))/(len(seq)+0.0))+1
    if (1 - stats.chi2.cdf(statistic, 2)) < significance:
        print("Statistic = "+str(statistic)+", Signficance = "+str((1 - stats.chi2.cdf(statistic, 2)))
              +", likely false.")
    else:
        print("Statistic = "+str(statistic)+", Signficance = "+str((1 - stats.chi2.cdf(statistic, 2)))
              +", likely true.")

#length of sequence must be >9
#chi square deg 1, non random is higher (one sided test)
def monobit(seq):
    global significance 
    n_0 = 0
    n_1 = 0
    n = len(seq)
    for bit in seq:
        if bit == 0:
            n_0 += 1
        else:
            n_1 += 1
    statistic = ((n_0-n_1)**2)/float(n)
    if (1 - stats.chi2.cdf(statistic, 1)) < significance:
        print("Statistic = "+str(statistic)+", Signficance = "+str((1 - stats.chi2.cdf(statistic, 1)))
              +", likely false.")
    else:
        print("Statistic = "+str(statistic)+", Signficance = "+str((1 - stats.chi2.cdf(statistic, 1)))
              +", likely true.")

def poker(seq,m):
    global significance
    k = math.floor(len(seq)/(m+0.0))
    if k < 5*(2**m):
        return "error, invalid m parameter."
    counter = [0 for i in range(2**m)]
    comb = list(itertools.product([0, 1], repeat=m))
    for i in range(k):
        tup = ()
        for j in range(m):
            tup = tup+(seq[(m*i)+j],)
        #tup = (seq[i],seq[(i)+1],seq[(i)+2],seq[(i)+3])
        for j in range(2**m):
            if (tup == comb[j]):
                counter[j] += 1
                break
    sum_counter = 0
    for i in range(2**m):
        sum_counter += (counter[i]**2)
    statistic = (((2**m)/(k+0.0))*(sum_counter))-k
    if (1 - stats.chi2.cdf(statistic, (2**m)-1)) < significance:
        print("Statistic = "+str(statistic)+", Signficance = "+str((1 - stats.chi2.cdf(statistic, (2**m)-1)))
              +", likely false.")
    else:
        print("Statistic = "+str(statistic)+", Signficance = "+str((1 - stats.chi2.cdf(statistic, (2**m)-1)))
              +", likely true.")

def runs(seq):
    n = len(seq)
    gap = 0
    block = 0
    zero = True
    k = 1
    while(True):
        e = (n-k+1+3)/(2**(k+1+2))
        if e >= 5:
            k += 1
        else:
            break
    e = (n-k+3)/(2**(k+2))  
    gap_counts = [0 for i in range(n)]
    block_counts = [0 for i in range(n)]
    for bit in seq:
        if ((zero == True) and (bit == 1)):
            zero = False
            if gap > 0:
                gap_counts[gap] += 1
                gap = 0
        if ((zero == False) and (bit == 0)):
            zero = True
            if block > 0:
                block_counts[block] += 1
                block = 0
        if bit == 0:
            gap += 1
        else:
            block += 1
            zero = False
    if gap > 0:
        gap_counts[gap] += 1
    if block > 0:
        block_counts[block] += 1
    sum_blocks = 0
    sum_gaps = 0
    for i in range(1,k+1):
        e_i = (n-i+3)/(2**(i+2))
        sum_blocks += ((block_counts[i] - e_i)**2)/e_i
        sum_gaps+= ((gap_counts[i] - e_i)**2)/e_i

    statistic = sum_blocks + sum_gaps
    if (1 - stats.chi2.cdf(statistic, (2*k)-2)) < significance:
        print("Statistic = "+str(statistic)+", Signficance = "+str((1 - stats.chi2.cdf(statistic, (2*k)-2)))
              +", likely false.")
    else:
        print("Statistic = "+str(statistic)+", Signficance = "+str((1 - stats.chi2.cdf(statistic, (2*k)-2)))
              +", likely true.")

def autocor(seq,d):
    if d > math.floor(len(seq)/2):
        return "error, sequence too small."
    xor_sum = 0
    for i in range(len(seq)-d):
        if (not(seq[i] and seq[i+d])) and (seq[i] or seq[i+d]):
            xor_sum +=1
    statistic = 2*(xor_sum-((len(seq)-d)/2.0))/math.sqrt(len(seq)-d)
    if (1 - stats.norm.cdf(statistic) < significance) or (stats.norm.cdf(statistic) < significance):
        print("Statistic = "+str(statistic)+", likely false.")
    else:
        print("Statistic = "+str(statistic)+", likely true.")


def universal(seq,L,Q,K):#Q should be at least 10*(2**L)
    #K hsould be at least 1000*(2**L)
    #sequence length should be at least 1010*(2**l)
    
    distribution_values = [
    [6, 5.2177052, 2.954],
    [7, 6.1962507 ,3.125],
    [8, 7.1836656 ,3.238],
    [9, 8.1764248 ,3.311],
    [10, 9.1723243 ,3.356],
    [11, 10.170032 ,3.384],
    [12, 11.168765 ,3.401],
    [13, 12.168070 ,3.410],
    [14, 13.167693 ,3.416],
    [15, 14.167488, 3.419],
    [16, 15.167379, 3.421]]

    if (not (L < 17) and (L > 5)) or (K<(1000*(2**L))) or (Q<10*(2**L)) or (len(seq) < L*1010*(2**L)) or (((K+Q)*L) > len(seq)):
        print(K,(1000*(2**L)))
        print(Q,10*(2**L))
        print(len(seq),(K+Q)*L)
        print("Parameter error")
    else:
        T = [0 for i in range(2**L)]
        for i in range(0,Q):
            T[to_int(seq[i*L:(i*L)+L])] = i
        total = 0
        for i in range(Q,Q+K):
            total += math.log(i-T[to_int(seq[i*L:(i*L)+L])],2)
            T[to_int(seq[i*L:(i*L)+L])] = i

        value = total/float(K)
        c = 0.7 - (0.8/float(L)) + ((4+(32/float(L)))*K**(-3/float(L)))/15.0
        stdev = c*math.sqrt(distribution_values[L-6][2]/float(K))
        statistic = abs((value - distribution_values[L-6][1])/(math.sqrt(2)*stdev))
        p_value = math.erfc(statistic)
        if (p_value < significance2):
            print("Statistic = "+str(statistic)+", Signficance = "+str(p_value)+", likely false.")
        else:
            print("Statistic = "+str(statistic)+", Signficance = "+str(p_value)+", likely true.")