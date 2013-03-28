#
# From http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank2fasta/
# Modified by Matthew Monroe on 9/19/2012
#


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

gbk_filename = "SourceFile.gbk"
faa_filename = "SourceFile_converted.fasta"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

processedCount = 0
for seq_record in SeqIO.parse(input_handle, "genbank") :
#    print "Dealing with GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="source" :
            record_organism = seq_feature.qualifiers['organism'][0]
        if seq_feature.type=="CDS" :
            processedCount = processedCount + 1
            if processedCount % 1000 == 0 :
                print "CDS entry %d: %s from %s" % (processedCount, seq_feature.qualifiers['locus_tag'][0], seq_record.id)

            assert len(seq_feature.qualifiers['translation'])==1

            # Note: Could use this to write out a standard protein fasta file entry
            # currentProteinRecord = SeqRecord(Seq(seq_feature.qualifiers['translation'][0], IUPAC.protein),
            #        id=seq_feature.qualifiers['locus_tag'][0], name="UnknownGene",
            #        description=seq_feature.qualifiers['product'][0])
            # SeqIO.write(currentProteinRecord, output_handle, "fasta")


            # Instead, we use this code so that we can include the organism name
            # Note: could include Locus using seq_record.name
            output_handle.write(">%s %s [%s]\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_feature.qualifiers['product'][0],
                   record_organism))

            residues = seq_feature.qualifiers['translation'][0]
            
            i = 0
            chunkSize = 60
            while (i < len(residues)):
                output_handle.write("%s\n" % (
                       residues[i:i+chunkSize]))
                i += chunkSize

output_handle.close()
input_handle.close()
