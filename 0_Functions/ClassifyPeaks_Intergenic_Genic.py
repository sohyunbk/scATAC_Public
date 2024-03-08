from pybedtools import BedTool
import pybedtools
import argparse

#conda activate  r_env
#~/.conda/envs/r_env/bin/python
###

def Make_Fai_Dic(GenomeSizeFile):
    infile = open(GenomeSizeFile,"r")
    Dic = {}
    for sLine in infile:
        if sLine.startswith("chr"):
            sList = sLine.strip().split("\t")
            Dic[sList[0]] = int(sList[1])
    return Dic
def ReadBedFiles_onlyChr(InfileName):
    infile = open(InfileName,"r")
    List=[]
    for sLine in infile:
        sList = sLine.strip().split("\t")
        if sList[0] in GenomeSize_Dic.keys() and int(sList[2]) < GenomeSize_Dic[sList[0]]:
            List.append(sLine)
    infile.close()
    BedtoolRead = BedTool(List)
    return BedtoolRead
def Change_PeakToSummit(PeakFileBedtools):
    NewList = []
    for i in open(PeakFileBedtools.fn):
        sList = i.strip().split("\t")
        nNewStart = int(sList[1])+250
        nNewEnd = int(sList[1])+251
        NewList.append(sList[0]+"\t"+str(nNewStart)+"\t"+str(nNewEnd)+"\t"+"_".join(sList)+"\n")
    BedtoolRead = BedTool(NewList)
    return BedtoolRead
def Outfile_forthe4thColumn(BedToolsInput,OutfileName):
    outfile = open(OutfileName,"w")
    for sLine in open(BedToolsInput.fn):
        OutLine = "\t".join(sLine.strip().split("\t")[3].split("_"))+"\n"
        outfile.write(OutLine)
    outfile.close()
def get_parser():
    parser = argparse.ArgumentParser(
        description="Call Peaks for scATAC data. \
    Requires cluster annnotations, as well as BED file ipput."
    )
    parser.add_argument(
        "-PeakFile",
        "--PeakFile",
        help="PeakFile",
        required=True,
        dest="peak",
    )
    parser.add_argument(
        "-BlackListFile",
        "--BlackListFile",
        help="BlackListFile",
        required=True,
        dest="BlackList",
    )

    parser.add_argument(
        "-FaiFile", "--FaiFile", help="GenomeSizefile", required=True, dest="Fai"
    )
    parser.add_argument(
        "-Ann_bed", "--Ann_bed", help="Ann_bed", required=True, dest="Ann_bed"
    )
    parser.add_argument(
        "-OutfilePath", "--OutfilePath", help="OutfilePath", required=True, dest="OutfilePath"
    )
    parser.add_argument(
        "-OutfileName", "--OutfileName", help="OutfileName", required=True, dest="OutfileName"
    )
    args = vars(parser.parse_args())
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    #CommonPeakFile =  "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/ComA619Bif3.unique500bpPeaks_sorted.bed"
    #BlackListFile = "/scratch/sb14489/0.Reference/Maize_B73/Zm.final_blaclist.Mito_Chloro_Chip.txt"
    #GenomeSize =  "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai"
    #Ann_bed = "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed"
    CommonPeakFile = args.peak
    BlackListFile = args.BlackList
    GenomeSize = args.Fai
    Ann_bed = args.Ann_bed
    GenomeSize_Dic = Make_Fai_Dic(GenomeSize)
    CommonPeak = ReadBedFiles_onlyChr(CommonPeakFile)
    BlackList = ReadBedFiles_onlyChr(BlackListFile)    # [1]
    Ann = ReadBedFiles_onlyChr(Ann_bed)
        #a.intersect(b, wo=True)
    CommonPeak_RemoveBlackList = CommonPeak-BlackList
    CommonSummit = Change_PeakToSummit(CommonPeak_RemoveBlackList)
    CommonPeak_Intergenic = CommonSummit-Ann
    CommonPeak_genic = CommonSummit + Ann
    Outfile_forthe4thColumn(CommonPeak_Intergenic,args.OutfilePath+"/"+args.OutfileName+"_Intergenic.bed")
    Outfile_forthe4thColumn(CommonPeak_genic,args.OutfilePath+"/"+args.OutfileName+"_Genic.bed")

    print("All")
    print(len(CommonPeak))
    print("After BlackList")
    print(len(CommonPeak_RemoveBlackList))
    print(len(CommonPeak_Intergenic) + len(CommonPeak_genic))
    print("Intergenic")
    print(len(CommonPeak_Intergenic))

#Intersect_Peak_Ann = CommonPeak_RemoveBlackList.intersect(Ann,wa=True, wb=True)
#Intersect_Peak_Ann = CommonPeak_RemoveBlackList.intersect(Ann,wo=True)

#Intersect_Peak_Ann_Merged = Intersect_Peak_Ann.merge(c="5,6,7,10", o="collapse,collapse,collapse,collapse")
#print(Intersect_Peak_Ann_Merged[30839])
#chr5	212650861	212651362	chr5,chr5	212650315,212647007	212651421,212651324	+,-
## Some ACRs overlapped with multiple genes. In this case, if they are ovelapped with at least one promoter, it is classified as Promoter-overlap
