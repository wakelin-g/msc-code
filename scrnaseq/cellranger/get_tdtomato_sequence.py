from Bio import GenBank

# genbank file from https://www.addgene.org/54642/sequences/
with open("./addgene-plasmid-54642-sequence-115699.gbk") as f:
    gbk = GenBank.read(f)

for i in gbk.features:
    for j in i.qualifiers:
        if "tdTomato" in j.value:
            start = int(i.location.split("..")[0])
            end = int(i.location.split("..")[1])

sequence = gbk.sequence[start - 1 : end]

header = ">tdTomato\n"

with open("./tdTomato.fa", "w") as f:
    f.write(header)
    f.write(sequence + "\n")

seqlen = len(sequence)

with open("./tdTomato.gtf", "w") as f:
    f.write("tdTomato\t")
    f.write("unknown\t")
    f.write("exon\t")
    f.write("1\t")
    f.write(f"{seqlen}\t")
    f.write(".\t")
    f.write("+\t")
    f.write(".\t")
    f.write('gene_id "tdTomato"; ')
    f.write('transcript_id "tdTomato"; ')
    f.write('gene_name "tdTomato"; ')
    f.write('gene_biotype "protein_coding"; ')
    f.write("\n")
