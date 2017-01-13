#####
## This script will find the frame with the fewest stop codons
#####
from Bio import SeqIO
from Bio.Alphabet import IUPAC, ProteinAlphabet
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord

src_filename = "SourceFile.fasta"
faa_filename = "SourceFile_converted.fasta"
src_filename = "BP_Sediment_Genomes_Jansson_Excerpt.fasta"
faa_filename = "BP_Sediment_Genomes_Jansson_Excerpt_converted_KeepBestProtein.fasta"
input_handle  = open(src_filename, "r")
output_handle = open(faa_filename, "w")

codonTableID = 1

def find_largest_polypeptide_in_DNA(seq, directionsToConsider="forward", tranlationTable=1):
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
            framePossibilities = [i[i.find("M"):] for i in trans.split("*") if "M" in i]
            allPossibilities += framePossibilities
            
    if directionsToConsider in ("reverse","both"):            
        # consider reverse complement DNA sequence as well
        # start translation from 1, 2 and 3 nucleotide
        for frame in range(3):
            trans = str(seq.reverse_complement()[frame:].translate(tranlationTable))
            framePossibilities = [i[i.find("M"):] for i in trans.split("*") if "M" in i]
            allPossibilities += framePossibilities

    # Find the length of each possible translated ORF
    allPossibilitiesLengths = [len(i) for i in allPossibilities]

    if len(allPossibilitiesLengths) == 0:
        raise Exception("no candidate ORFs")

    # Select the longest translated ORF
    proteinAsString = allPossibilities[allPossibilitiesLengths.index(max(allPossibilitiesLengths))]

    return Seq(proteinAsString, alphabet=ProteinAlphabet)

for seq_record in SeqIO.parse(input_handle, "fasta", alphabet=IUPAC.ambiguous_dna):
    try:
        largestProteinSeq = find_largest_polypeptide_in_DNA(seq_record.seq,"forward", codonTableID)
        largestProteinRecord = SeqRecord(largestProteinSeq, seq_record.name)
        largestProteinRecord.description = seq_record.description

        SeqIO.write(largestProteinRecord, output_handle, "fasta")

    except Exception as inst:
        print("Error translating {0}, {1}".format(seq_record.name, inst.args[0]))

input_handle.close()
output_handle.close()
