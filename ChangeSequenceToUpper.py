from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord

src_filename = "SourceFile_Lowercase.fasta"
faa_filename = "SourceFile_Uppercase.fasta"
input_handle  = open(src_filename, "r")
output_handle = open(faa_filename, "w")

i = 0
for seq_record in SeqIO.parse(input_handle, "fasta", alphabet=IUPAC.ambiguous_dna) :
    i = i + 1
    if i % 100 == 0 :
        print("Entry {0}: {1}".format(i, seq_record.name))
   
    try:
        proteinRecord = SeqRecord(seq_record.seq.upper(), seq_record.name)
        proteinRecord.description = seq_record.description
        
        SeqIO.write(proteinRecord, output_handle, "fasta")
       
    except Exception as inst:
        print("Error changing case for {0}, {1}".format(seq_record.name, inst.args[0]))
        output_handle.write(">%s %s\n%s\n" % (
           seq_record.name,
           seq_record.description,
           seq_record.seq))
        pass

output_handle.close()
input_handle.close()
