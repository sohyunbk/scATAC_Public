import os,sys, glob

WD="/scratch/sb14489/3.scATAC_flo/2.Mapped/"

#List=["1_A619_2","2_rel2_2","3_bif3_2","4_relk1_2"]
List=["1_A619_2","2_rel2_2"]
Outfile = open(WD+"MappedSummary_2ndReplicates.csv","w")
k = 0

 
for i in List:
    File = WD+i+"/outs/summary.csv"
    infile = open(File,"r")
    FirstLine = infile.readline()
    SecondLine = infile.readline()
    if k == 0:
        Outfile.write("SampleID,"+FirstLine)
        Outfile.write(i+","+SecondLine)
    else:
        Outfile.write(i+","+SecondLine)
    k=+1
    infile.close()

Outfile.close()
