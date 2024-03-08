import os,sys, glob,gzip

## This one should be done after BarcodeFixing

List=["1_A619","1_A619_2","2_rel2","2_rel2_2","3_bif3","3_bif3_2","4_relk1","4_relk1_2"]
#List=["1_A619_2","2_rel2_2"]
#List=["axillarybud1","axillarybud2","crownRoot1","crownRoot2"]
#List=["Tassel_ba1", "Tassel_bif1", "Seedling_IAA2"]
OutputName="MappingSummary_Ear"

################################################################################
################################################################################

WD="/scratch/sb14489/3.scATAC_flo/"

Outfile = open(WD+"/2.Mapped/"+OutputName+".csv","w")

k = 0 
Outfile.write("SampleName,Total#ofReads,RatioandNumberOfMappedReads,UniqueBarcodeNumber\n")

for i in List:
    infile = open(WD+"/2.Mapped/"+i+"_MappedReads.txt","r")
    for sLine in infile:
        #sList = sLine.strip().split("\t")
        if "paired in sequencing" in  sLine:
            Total = int(sLine.split(" ")[0])/2
        if "properly paired" in sLine:
            Mapped = int(sLine.split(" ")[0])/2
            MappingRate = sLine.split("(")[1].split(" ")[0]
    infile.close()
    infile2 = open(WD+"/4.Bam_FixingBarcode/"+i+"_Unique.bed","r")

    Dic = {}
    for sLine in infile2:
        sList = sLine.strip().split("\t")
        Barcode = sList[3]
        Dic.setdefault(Barcode,"")

    infile2.close()
    #print(sFiles)
    nUnique = len(Dic.keys())
    infile2.close()
    Outfile.write(i+","+str(Total)+","+str(Mapped)+" ("+str(MappingRate)+"),"+str(nUnique)+"\n")


Outfile.close()
