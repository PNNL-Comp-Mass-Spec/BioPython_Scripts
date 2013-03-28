#####
## This script will find the longest series of amino acids that starts with "M" and has no stop codons 
#####
from Bio import SeqIO
from Bio.Alphabet import IUPAC, ProteinAlphabet
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord

src_filename = "SourceFile.fasta"
faa_filename = "SourceFile_converted.fasta"
src_filename = "BP_Sediment_Genomes_Jansson_Excerpt.fasta"
faa_filename = "BP_Sediment_Genomes_Jansson_Excerpt_converted_KeepBestFrame.fasta"
input_handle  = open(src_filename, "r")
output_handle = open(faa_filename, "w")

codonTableID = 1

def find_best_frame_in_DNA(seq, directionsToConsider="forward", tranlationTable=1):
    """
    directionToConsider:
        forward - normal DNA direction
        reverse - reverse complement
        both - both of the abover
    """
    allPossibilities = []
    if directionsToConsider in ("forward","both"):
        # start translation from 1, 2 and 3 nucleotide
        for frame in range(3):
            trans = str(seq[frame:].translate(tranlationTable))
            allPossibilities.append(trans)
                        
    if directionsToConsider in ("reverse","both"):            
        # consider reverse complement DNA sequence as well
        # start translation from 1, 2 and 3 nucleotide
        for frame in range(3):
            trans = str(seq.reverse_complement()[frame:].translate(tranlationTable))
            allPossibilities.append(trans)

    # Count the number of stop codons in each frame
    stopCodonCount = [i.count("*") for i in allPossibilities]

    # Select the frame with the fewest stop codons
    bestFrame = allPossibilities[stopCodonCount.index(min(stopCodonCount))]
    
    return Seq(bestFrame, alphabet=ProteinAlphabet)

for seq_record in SeqIO.parse(input_handle, "fasta", alphabet=IUPAC.ambiguous_dna):
    try:
        bestFrame = find_best_frame_in_DNA(seq_record.seq,"forward", codonTableID)
        bestFrameRecord = SeqRecord(bestFrame, seq_record.name)
        bestFrameRecord.description = seq_record.description

        SeqIO.write(bestFrameRecord, output_handle, "fasta")

    except Exception as inst:
        print "Error translating %s, %s" % (seq_record.name, inst.args[0])

input_handle.close()
output_handle.close()
