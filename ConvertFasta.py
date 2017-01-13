from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord

src_filename = "SourceFile.fasta"
faa_filename = "SourceFile_converted.fasta"
src_filename = "BP_Sediment_Genomes_Jansson_Excerpt.fasta"
faa_filename = "BP_Sediment_Genomes_Jansson_Excerpt_converted.fasta"
input_handle  = open(src_filename, "r")
output_handle = open(faa_filename, "w")

codonTableID = 1

# Could use this to load all records into memory
# records = list(SeqIO.parse(src_filename, "fasta"))

i = 0
for seq_record in SeqIO.parse(input_handle, "fasta", alphabet=IUPAC.ambiguous_dna) :
    i = i + 1

    if i % 100 == 0 :
        print("Entry {0}: {1}".format(i, seq_record.name))
    
    # Use this to see help on the translate method:
    # help (records[0].seq.translate)
    # You could specify the translation table like this:
    # seq_record.seq.translate(table="Bacterial")
    # Note that using cds="true" instructs the code to verify that each sequence is a true CDS
   
    try:
        proteinRecord = SeqRecord(seq_record.seq.translate(cds="false", to_stop="false", table=codonTableID), seq_record.name)
        proteinRecord.description = seq_record.description
        
        SeqIO.write(proteinRecord, output_handle, "fasta")
                
        #output_handle.write(">%s %s\n%s\n" % (
        #   seq_record.name,
        #   seq_record.description,
        #   seq_record.seq.translate(cds="false", to_stop="false", table=codonTableID)))
           
    except Exception as inst:
        print("Error translating {0}, {1}".format(seq_record.name, inst.args[0]))

        proteinRecord = SeqRecord(translate(seq_record.seq, codonTableID), seq_record.name)
        proteinRecord.description = seq_record.description + " (translation warning: " + inst.args[0] + ")"
        
        SeqIO.write(proteinRecord, output_handle, "fasta")
        
        pass

output_handle.close()
input_handle.close()