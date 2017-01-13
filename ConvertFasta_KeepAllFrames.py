#####
## This script will translate all three or all six frames
## For three frames, set translate_direction to "forward" or "reverse"
## For six frames,   set translate_direction to "both"
#####
from Bio import SeqIO
from Bio.Alphabet import IUPAC, ProteinAlphabet
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord

src_filename = "SourceFile.fasta"
faa_filename = "SourceFile_converted.fasta"
translate_direction = "both"
minimum_peptide_length = 10

src_filename = "BP_Sediment_Genomes_Jansson.fasta"
faa_filename = "BP_Sediment_Genomes_Jansson_converted_KeepAllFrames.fasta"

input_handle  = open(src_filename, "r")
output_handle = open(faa_filename, "w")

#codonTableID = 1      # Use for eukaryotes
codonTableID = 11     # Use for bacteria


def translate_and_write_DNA_frames(seq, directionsToConsider="forward", tranlationTable=1):
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

    i = 0
    for currentFrame in allPossibilities:
        i = i + 1
        currentProtein = Seq(currentFrame, alphabet=ProteinAlphabet)

        if len(currentProtein) >= minimum_peptide_length:
            currentProteinRecord = SeqRecord(currentProtein, seq_record.name)
            currentProteinRecord.id = currentProteinRecord.id + "." + str(i)
            currentProteinRecord.description = seq_record.description + "; frame " + str(i)
    
            SeqIO.write(currentProteinRecord, output_handle, "fasta")
        

    return
    

for seq_record in SeqIO.parse(input_handle, "fasta", alphabet=IUPAC.ambiguous_dna):
    try:
        translate_and_write_DNA_frames(seq_record.seq, translate_direction, codonTableID)
    
    except Exception as inst:
        print("Error translating {0}, {1}".format(seq_record.name, inst.args[0]))

input_handle.close()
output_handle.close()
