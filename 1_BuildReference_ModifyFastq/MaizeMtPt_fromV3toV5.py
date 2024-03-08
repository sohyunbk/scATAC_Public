WD="/scratch/sb14489/3.scATAC/1.MaizeExample/0.Reference/"

## 220901 -- Updated --> Remove all the scaffolds
############ Fa file ###################################

Fa_v4 = open(WD+"v3/B73_RefGen_v3.fa","r") ## it's V3

Dic = {}
for sLine in Fa_v4:
    if sLine.startswith(">"):
        Key = sLine
        Dic.setdefault(Key,[])
    else:
        Dic[Key].append(sLine)

Fa_v4.close()

Fa_v5 = open(WD+"v5/Zm-B73-REFERENCE-NAM-5.0.fa","r")
Fa_Write = open(WD+"/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd.fa","w")
for sLine in Fa_v5:
    Fa_Write.write(sLine)

Fa_v5.close()
print(Dic.keys())

#>Mt dna:chromosome chromosome:AGPv3:Mt:1:569630:1
#>Pt dna:chromosome chromosome:AGPv3:Pt:1:140384:1

Fa_Write.write(">Mt\n")
for sLine in Dic[">Mt dna:chromosome chromosome:AGPv3:Mt:1:569630:1\n"]:
    Fa_Write.write(sLine)
Fa_Write.write(">Pt\n")
for sLine in Dic[">Pt dna:chromosome chromosome:AGPv3:Pt:1:140384:1\n"]:
    Fa_Write.write(sLine)
Fa_Write.close()

############ Gff file ##################################


Fa_v4 = open(WD+"v3/Zea_mays.AGPv3.22.gff3","r") ## it's V3
Dic = {}
for sLine in Fa_v4:
    if sLine.startswith("chrMt") or sLine.startswith("chrPt"):
        Key = sLine.strip().split("\t")[0]
        Dic.setdefault(Key,[])
        Dic[Key].append(Key.replace("chr","")+"\t"+"\t".join(sLine.split("\t")[1:len(sLine.split("\t"))]))

Fa_v4.close()

Fa_v5 = open(WD+"v5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3","r")
Fa_Write = open(WD+"Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd.gff3","w")
for sLine in Fa_v5:
    Fa_Write.write(sLine)

Fa_v5.close()
print(Dic.keys())

#>Mt dna:chromosome chromosome:AGPv3:Mt:1:569630:1
#>Pt dna:chromosome chromosome:AGPv3:Pt:1:140384:1

for i in Dic.keys():
    print(i)
    for k in Dic[i]:
        Fa_Write.write(k)
Fa_Write.close()

#-rw-r--r--. 1 sb14489 rjslab  159980738 Aug 15 11:14 Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd.gff3
#-rw-r--r--. 1 sb14489 rjslab 2210080866 Aug 15 11:14 Zm-B73-REFERENCE-NAM-5.0_MtPtAdd.fa
