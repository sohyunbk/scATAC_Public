cd /scratch/sb14489/3.scATAC/1.MaizeExample/0.Reference/v5
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
gzip -d Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
gzip -d Zm-B73-REFERENCE-NAM-5.0.fa.gz

cd /scratch/sb14489/3.scATAC/1.MaizeExample/0.Reference/v4
wget https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz
gzip -d Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz
gzip -d Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz

cd /scratch/sb14489/3.scATAC/1.MaizeExample/0.Reference/v3
wget https://download.maizegdb.org/B73_RefGen_v3/B73_RefGen_v3.fa.gz
wget https://download.maizegdb.org/B73_RefGen_v3/Zea_mays.AGPv3.22.gff3.gz
gzip -d B73_RefGen_v3.fa.gz
gzip -d Zea_mays.AGPv3.22.gff3.gz
