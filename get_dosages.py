import sys
import glob

# sys.argv[1] = directory to dosages containing files *.matrix.maf5 split by chromosome.
# sys.argv[2] = file which contains col 1 = SNP, col2 = gene
# sys.argv[3] = file out to write SNP dosages

dosage_files=glob.glob(sys.argv[1]+'*.matrix.maf5')

eqtls={}
with open(sys.argv[2],'r') as f:
   for line in (line.strip().split() for line in f):
      eqtls[line[1]]=line[0]

snps=set(eqtls.values())

with open(sys.argv[3],'w') as output:
   for x in dosage_files:
      with open(x,'r') as results:
         for line in (line.strip().split() for line in results):
            if line[0] in snps:
               output.write("\t".join(str(x) for x in line)+"\n")
