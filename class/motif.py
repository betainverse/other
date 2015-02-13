#!/usr/bin/python
"""
motif.py
Given: DNA sequence

Return: Positions of 10 motifs
Given: non-motif sequences have C->G transition of 0.08
motif sequences have A->T transition of 0.7 and T->G of 0.9
"""
from operator import mul
import pylab

# Transition probabilities, ignoring begin and end

Tr_bg = {}
Tr_bg['AA']=Tr_bg['AC']=Tr_bg['AG']=Tr_bg['AT']=0.25
Tr_bg['GA']=Tr_bg['GC']=Tr_bg['GG']=Tr_bg['GT']=0.25
Tr_bg['TA']=Tr_bg['TC']=Tr_bg['TG']=Tr_bg['TT']=0.25
Tr_bg['CG']=0.08
Tr_bg['CA']=Tr_bg['CC']=Tr_bg['CT']=(1-Tr_bg['CG'])/3.0

Tr_mot = {}
Tr_mot['AT']=0.7
Tr_mot['AA']=Tr_mot['AC']=Tr_mot['AG']=(1-Tr_mot['AT'])/3.0
Tr_mot['CA']=Tr_mot['CC']=Tr_mot['CG']=Tr_mot['CT']=0.25
Tr_mot['GA']=Tr_mot['GC']=Tr_mot['GG']=Tr_mot['GT']=0.25
Tr_mot['TG']=0.9
Tr_mot['TA']=Tr_mot['TC']=Tr_mot['TT']=(1-Tr_mot['TG'])/3.0

# Given a list of numbers, return the product of those numbers
def list_product(list_of_num):
    return reduce(mul,list_of_num,1)

# Given a sequence and a dictionary of transition probabilities,
# calculate the probability of producing that sequence using that dictionary
def calc_prob(seq,tr):
    transition_probs = [tr[seq[i:i+2]] for i in range(len(seq)-2)]  
    return list_product(transition_probs)

# Given a sequence, calculate the likelihood ratio for motif/background
def motif_likelihood_ratio(seq):
    global Tr_mot
    global Tr_bg
    return calc_prob(seq,Tr_mot)/calc_prob(seq,Tr_bg)

## # Given a sequence, calculate the likelihood ratio that each subsequence
## # of length n is a motif or not; return the 10 most likely positions. 
## # This dictionary seems backwards.
## # Recall that positions start at 0.
## def calc_subsequence_motif(seq,n):
##     best_ratios = {}
##     for i in range(len(seq)-n+1):
##         ratio = motif_likelihood_ratio(seq[i:i+n])
##         if len(best_ratios.values()) < 10:
##             best_ratios[ratio]=i
##         elif ratio > max(best_ratios.keys()):
##             del best_ratios[min(best_ratios.keys())]
##             best_ratios[ratio]=i
##         print seq[i:i+n], ratio
##     most_likely_motifs={}
##     for position in best_ratios.values():
##         most_likely_motifs[position] = seq[i:i+n]
##     return most_likely_motifs

# Given a sequence, calculate the likelihood ratio that each subsequence
# of length n is a motif or not. Recall that positions start at 0.
# Return a dictionary of likelihood ratios and sequences for each position
def calc_subsequence_motif(seq,n):
    ratios = {}
    for i in range(len(seq)-n+1):
        subsequence = seq[i:i+n]
        ratios[i] = motif_likelihood_ratio(subsequence)
    return ratios

# Given a dictionary of likelihood ratios, select
# the subsequence with the highest likelihood of being a motif. Return
# the index of that subsequence, along with the likelihood and the sub-
# sequence itself.

def get_best_motif(ratio_dict):
    max_likelihood = 0
    best_index = 0
    for key in ratio_dict.keys():
        if ratio_dict[key]>max_likelihood:
            max_likelihood = ratio_dict[key]
            best_index = key
    return key,max_likelihood

# Determine whether two motifs overlap, given their indices and length.
# Return true if they are separate, and do not overlap; false otherwise.
# Here is where we can adjust to force a gap between the motifs. 

def is_separate(index1,index2,length):
    if abs(index1-index2)<length:
        return False
    else:
        return True

# Determine whether a motif overlaps with any others in a list, given
# its index and a list of the indices of the others, and the length of motifs.

def is_unique(index1,indices,length):
    if len(indices)<1:
        return True
    for i in indices:
        if not is_separate(index1,i,length):
            return False
    return True

# Given a set of subsequences and their likelihood ratios, cull the list
# to obtain the 10 best candidate motifs. If two sequences in the top 10
# have significant overlap, drop the lower. For now, define significant
# overlap as any overlap at all. However, we currently will allow motifs
# to be adjacent without any intervening sequence, although this doesn't
# really make sense. 

def cull_motifs(ratio_dict,motif_length):
    top10 = {}
    for i in range(10):
        max_likelihood = 0
        best_index = 0
        for key in ratio_dict.keys():
            if ratio_dict[key]>max_likelihood:
                dist_from_others = [abs(key-akey) for akey in top10.keys()]
                unique = True
                for dist in dist_from_others:
                    if dist < motif_length:
                        unique = False
                if unique:
                    max_likelihood = ratio_dict[key]
                    best_index = key
        top10[best_index] = max_likelihood
    return top10

def print_results_a(most_likely_motifs,seq,length):
    for position in most_likely_motifs.keys():
        print position, most_likely_motifs[position],seq[position:position+length]

def get_sequence(filepath):
    openfile = open(filepath,'r')
    seq = openfile.readline()
    openfile.close()
    return seq

def plot_ratios(ratio_dict):
    xlist = ratio_dict.keys()
    ylist = [ratio_dict[key] for key in ratio_dict.keys()]
    fig=pylab.figure()
    pylab.scatter(xlist,ylist,marker='o',s=5)
    pylab.xlabel("sequence index")
    pylab.ylabel("odds ratio")
    axes = fig.gca()
    axes.set_yscale('log')
#    axes.set_aspect('equal')
    pylab.xlim(min(xlist),max(xlist))
    pylab.ylim(min(ylist),max(ylist))
    pylab.savefig('motifs.pdf')

def main():
    
    global Tr_bg
    global Tr_mot
    #print "list_product([1,2,3,4,5]) =",list_product([1,2,3,4,5])
    #print "Pr(AAAA) =", calc_prob('AAAA',Tr_bg)
    #subseq = 'ATGGCTGGTGTGCATGATGTGATGATGCTGAATGCCATGCGCCTGTGTGT'
    #print motif_likelihood_ratio('ATGGCTGGTGTGCATGATGTGATGATGCTGAATGCCATGCGCCTGTGTGT')
    #print "likelihood(AAAA) =", motif_likelihood_ratio('AAAA')
    #print "likelihood(ATAT) =", motif_likelihood_ratio('ATAT')
    #print "likelihood ratios for AATGCCTGTGCCGTAAAGTCCA:"
    #print calc_subsequence_motif('AATGCCTGTGCCGTAAAGTCCA',10)
    #print get_sequence('new_seq_1.txt')[-20:]
    seq = get_sequence('new_seq_1.txt')
    print_results_a(cull_motifs(calc_subsequence_motif(seq,50),50),seq,50)
    #print cull_motifs(calc_subsequence_motif(get_sequence('new_seq_1.txt'),50),50)
    #plot_ratios(calc_subsequence_motif(get_sequence('new_seq_1.txt'),50))
    
main()
