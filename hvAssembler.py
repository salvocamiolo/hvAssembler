import os,sys
from Bio import SeqIO
from Bio import Seq


read1 = sys.argv[1]
read2 = sys.argv[2]
condaDir = sys.argv[3]
outputFolder = sys.argv[4]

os.system("mkdir -p "+outputFolder)

kmerReadDict = {}

fq1 = {}
fq2 = {}
#Load reads kmers in memory
print("Loading kmer for read1 in memory....")
for seq_record in SeqIO.parse(read1,"fastq"):
    locus = (str(seq_record.id).split("/"))[0]
    sequence = str(seq_record.seq)
    if not locus in fq1:
        fq1[locus] = seq_record
    for a in range(0,len(sequence)-17,+17):
        kmer = sequence[a:a+17]
        if not kmer in kmerReadDict:
            kmerReadDict[kmer] = set()
        kmerReadDict[kmer].add(locus)

print("Loading kmer for read2 in memory....")
for seq_record in SeqIO.parse(read1,"fastq"):
    locus = (str(seq_record.id).split("/"))[0]
    sequence = str(seq_record.seq)
    if not locus in fq2:
        fq2[locus] = seq_record
    for a in range(0,len(sequence)-17,+17):
        kmer = sequence[a:a+17]
        if not kmer in kmerReadDict:
            kmerReadDict[kmer] = set()
        kmerReadDict[kmer].add(locus)


#Load specific kmers in memory
kmerDict = {}
infile = open("mainDB_seqs_filtered.txt")
infile.readline()
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if not fields[0] in kmerDict:
        kmerDict[fields[0]] = {}
    if not fields[1] in kmerDict[fields[0]]:
        kmerDict[fields[0]][fields[1]] = []
    kmerDict[fields[0]][fields[1]] = fields[2]

infile.close()

#Scan hypervariable gene kmers
outfile = open("./"+outputFolder+"/foundGenotypes.txt","w")
outfile.write("Gene\tGenotype\tNum_Found_Kmers\tNum_Found_reads")
genotyeSpecificReads = {}
numFoundKmers = {}
numFoundReads = {}
for gene in kmerDict:
    print("Retrieving reads for genotype %s" %gene)
    for genotype in kmerDict[gene]:
        if not genotype in numFoundKmers:
            numFoundKmers[genotype] = 0
        if not genotype in numFoundReads:
            numFoundReads[genotype] = 0
        for kmer in (kmerDict[gene][genotype]).split(","):
            if kmer in kmerReadDict:
                numFoundKmers[genotype]+=1
                numFoundReads[genotype]+=len(kmerReadDict[kmer])
                if not (gene,genotype) in genotyeSpecificReads:
                    genotyeSpecificReads[(gene,genotype)] = set()
                for read in kmerReadDict[kmer]:
                    genotyeSpecificReads[(gene,genotype)].add(read)
    

    for genotype in numFoundKmers:
        if numFoundKmers[genotype]>0:
            print("%d reads were found in genotype %s of gene %s" %(numFoundKmers[genotype],genotype,gene))
            outfile.write(gene+"\t"+genotype+"\t"+str(numFoundKmers[genotype])+"\t"+str(numFoundReads[genotype])+"\n")
        

os.system("mkdir -p ./"+outputFolder+"/kmerSpecificReads")
os.system("mkdir -p ./"+outputFolder+"/kmerSpecificScaffolds")
for item in genotyeSpecificReads:
    print("Assembling reads for gene/genotype %s/%s" %(item[0],item[1]))
    newRead1 = "./"+outputFolder+"/kmerSpecificReads/"+item[0]+"_"+item[1]+"_read_1.fastq"
    newRead2 = "./"+outputFolder+"/kmerSpecificReads/"+item[0]+"_"+item[1]+"_read_2.fastq"
    outfq1 = open(newRead1,"w")
    outfq2 = open(newRead2,"w")

    for readName in genotyeSpecificReads[item]:
        SeqIO.write(fq1[readName],outfq1,"fastq")
        SeqIO.write(fq2[readName],outfq2,"fastq")

    
    
    outfq1.close()
    outfq2.close()

    os.system(condaDir+"/bin/fastq_to_fasta -i "+newRead1+ " -o read1.fasta")
    os.system(condaDir+"/bin/fastq_to_fasta -i "+newRead2+ " -o read2.fasta")
    os.system("cat read1.fasta read2.fasta >all.fasta")
    

    os.system(condaDir+"/bin/cap3  all.fasta >out 2>&1")
    if os.path.isfile("all.fasta.cap.contigs") == True:
        os.system("cp all.fasta.cap.contigs ./"+outputFolder+"/kmerSpecificScaffolds/"+item[0]+"_"+item[1]+"_scaffolds.fasta")
    os.system("rm -rf all.fasta* read1.fasta read2.fasta null")


    





        




