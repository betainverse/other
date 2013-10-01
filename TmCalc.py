#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
"""
TmCalc.py -s <primer sequence> [-m <mismatch sequence>]
calculates the melting temperature of a primer, based on its length, GC
content, and number of mismatches, using the formula from the Stratagene
Quikchange kit: Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch

 * N is the primer length in bases
 * values for %GC and %mismatch are whole numbers

Primers should be between 25 and 45 bases in length, with a melting
temperature (Tm) of >= 78C. Primers longer than 45 bases may be used,
but using longer primers increases the likelihood of secondary structure
formation, which may affect the efficiency of the mutagenesis reaction.

This script should be modified to make better use of Biopython Seq.

"""

### Import some libraries #################################################

# OptionParser helps us parse the command line
from optparse import OptionParser
from Bio.Seq import Seq

### Helper Functions ######################################################

### Main body of the script ###############################################

def main():
    
    parser = OptionParser() # creates an instance of the parser
    parser.usage = "%prog -s primer_sequence [-m mismatch_sequence]"
    parser.description = "%prog calculates the melting temperature of a primer, based on the formula from the Stratagene Quikchange manual, Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch, given the sequence of the primer and the sequence of any mismatched regions of the primer."
    parser.epilog = "Example: TmCalc.py -s ggtGGAGTATATACCCTGTGAGAGAGGAGTCGGaagcc -m T"

    parser.add_option("-s", "--sequence", dest="sequence",type="string",default=None,
                      help="Primer sequence.")
    parser.add_option("-m", "--mismatch", dest="mismatch",type="string",default=None,
                      help="Sequence of any mismatched regions. Feel free to concatenate two separate mismatching regions. Be sure to accurately represent the total number of mismatched bases, and the total number of mismatched Gs and Cs.")
    parser.add_option("-n", "--name", dest="name",type="string",default=None,
                      help="You can provide a name for your primer, and it will appear with the results. A name is not required.")

    # Now parse the command-line options
    (options, args) = parser.parse_args()
    sequence = options.sequence
    mismatch = options.mismatch
    name = options.name

    if sequence == None:
        parser.print_help()
        parser.error("Please enter a primer sequence.")

    if mismatch == None:
        mismatch = ""

    if args:
        parser.error("If you have spaces in your sequence, you must enclose your sequence in quotes. Please try again. The following sequences were not parsed properly: %s"%args)

# Remove spaces in either sequence, convert to upper case
    subseqs = sequence.split()
    nwsSequence = ''.join(subseqs).upper()
    mismatchSubseqs = mismatch.split()
    nwsMismatch = ''.join(mismatchSubseqs).upper()

    # Replace U with T.

    ATGCsequence = nwsSequence.replace('U','T')
    ATGCmismatch = nwsMismatch.replace('U','T')

    # Count occurrences of A,T,G,C
    
    As = ATGCsequence.count('A')
    Ts = ATGCsequence.count('T')
    Gs = ATGCsequence.count('G')
    Cs = ATGCsequence.count('C')
    
    Total = As + Ts + Gs + Cs
    Length = len(ATGCsequence)

    mAs = ATGCmismatch.count('A')
    mTs = ATGCmismatch.count('T')
    mGs = ATGCmismatch.count('G')
    mCs = ATGCmismatch.count('C')

    mTotal = mAs + mTs + mGs + mCs
    mLength = len(ATGCmismatch)

    if Total != Length or mTotal != mLength:
        parser.error("Your sequence may contain only A,G,C,T, or U, in upper case or lower case, and white space. Any other characters are not permitted.")

    mGC = mGs + mCs
    GC = Gs + Cs - mGC
    
    percentGC = 100*GC/Total
    lengthUnder675 = 675.0/Total
    percentMismatch = 100*mTotal/Total
    Tm = 81.5 + (0.41*percentGC) - lengthUnder675 - percentMismatch 

    if name:
        print ""
        print name

    print ""
    print "Primer sequence: %s"%sequence
    print "Mismatch sequence: %s"%mismatch
    print "Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch"
    print "GC: %d"%GC, "%GC:", "%d"%percentGC
    print "N: %d\t675/N: %4.2f"%(Total,lengthUnder675)
    print "%mismatch:", "%d"%percentMismatch
    print "Tm = %5.2f"%Tm
    
    dna = Seq(ATGCsequence)
    print "Sequence: ", dna
    print "Reverse Complement: ", dna.reverse_complement()
    synthfw=str(dna).replace('A','5').replace('C','6').replace('G','7').replace('T','8')
    synthrev=str(dna.reverse_complement()).replace('A','5').replace('C','6').replace('G','7').replace('T','8')
    print "Synthesizer-fw: ",synthfw
    print "Synthesizer-rev: ",synthrev


# Execute everything
main()

