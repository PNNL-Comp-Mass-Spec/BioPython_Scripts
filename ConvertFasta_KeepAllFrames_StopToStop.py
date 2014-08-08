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
translate_direction = "forward"
minimum_peptide_length = 10

src_filename = "PTLD35.fasta"
faa_filename = "PTLD35_converted_stop-to-stop_fwd.fasta"

input_handle  = open(src_filename, "r")
output_handle = open(faa_filename, "w")

codonTableID = 1      # Use for eukaryotes
#codonTableID = 11     # Use for bacteria


def translate_and_write_DNA_frames(seq, directionsToConsider="forward", tranlationTable=1):
    """
    directionToConsider:
        forward - normal DNA direction
        reverse - reverse complement
        both - both of the abover
    """
    allPossibilities = []
    allPossibilitiesFrameInfo = []
    
    if directionsToConsider in ("forward","both"):
        # start translation from 1, 2 and 3 nucleotide
        currentFrameNum = 0
        for frame in range(3):
            currentFrameNum = currentFrameNum + 1
            trans = str(seq[frame:].translate(tranlationTable))
            currentPeptide = 0
            for i in trans.split("*"):
                currentPeptide = currentPeptide + 1
                allPossibilities.append(i)
                allPossibilitiesFrameInfo.append("fwd_" + str(currentFrameNum) + "_" + str(currentPeptide))
                        
    if directionsToConsider in ("reverse","both"):            
        # consider reverse complement DNA sequence as well
        # start translation from 1, 2 and 3 nucleotide
        currentFrameNum = 0
        for frame in range(3):
            currentFrameNum = currentFrameNum + 1
            trans = str(seq.reverse_complement()[frame:].translate(tranlationTable))
            currentPeptide = 0
            for i in trans.split("*"):
                currentPeptide = currentPeptide + 1
                allPossibilities.append(i)
                allPossibilitiesFrameInfo.append("rev_" + str(currentFrameNum) + "_" + str(currentPeptide))

    indexInfo = -1
    for currentFrame in allPossibilities:
        indexInfo = indexInfo + 1
        currentProtein = Seq(currentFrame, alphabet=ProteinAlphabet)

        if len(currentProtein) >= minimum_peptide_length:
            currentProteinRecord = SeqRecord(currentProtein, seq_record.name)
            currentProteinRecord.id = currentProteinRecord.id + "_" + allPossibilitiesFrameInfo[indexInfo]
            currentProteinRecord.description = seq_record.description + "; frame " + allPossibilitiesFrameInfo[indexInfo][4] + " " + allPossibilitiesFrameInfo[indexInfo][0:3]
    
            SeqIO.write(currentProteinRecord, output_handle, "fasta")
        

    return
    

for seq_record in SeqIO.parse(input_handle, "fasta", alphabet=IUPAC.ambiguous_dna):
    try:
        translate_and_write_DNA_frames(seq_record.seq, translate_direction, codonTableID)
    
    except Exception as inst:
        print("Error translating %s, %s" % (seq_record.name, inst.args[0]))

input_handle.close()
output_handle.close()
