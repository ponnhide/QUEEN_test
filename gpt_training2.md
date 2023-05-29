# Practial example codes of QUEEN 

By using QUEEN package, you can simulate any DNA cloning process and construct the DNA sequence of a taget DNA construct. Here, I will introduce some practical example codes to produce generala operation processes in DNA cloning experiments.  

## Example 1: pRS112 construction  
QUEEN script for lentiviral GFP expression plasmid construction by RE digestion and ligaiton process.
1. An enhanced green fluorescent protein (eGFP)-encoding cassette was amplified from the pLV-eGFP plasmid using the primer pair RS204/SI627 that has overhang sequences encoding EcoRI and BamHI restriction digestion sites. 
2. The eGFP-encoding cassette was digested by EcoRI and BamHI.
3. The destination plasmid pLVSIN-CMV-Puro (Takara) was also digested by EcoRI and BamHI.
4. The eGFP-encoding cassette was cloned into the destination vector using T4 DNA ligase (NEB).

```Python
from QUEEN.queen import *
from QUEEN import cutsite as cs
#Load pLVSIN_CMV_pur plasmid object from a Benchling share link
QUEEN(record="https://benchling.com/s/seq-0r6kj1kOMYhArIFPTIg5", dbtype="benchling", product="pLVSIN_CMV_pur") 

#Load pLV_eGFP plasmid object using an Addgene ID
QUEEN(record="36083", dbtype="addgene", product="pLV_eGFP")

pn = "PCR" 
pd = "1. An enhanced green fluorescent protein (eGFP)-encoding cassette was amplified from the pLV-eGFP plasmid using the primer pair RS204/SI627 that has overhang sequences encoding EcoRI and BamHI restriction digestion sites."
RS204 = QUEEN(seq="TCCGGTGAATTCCCGAGCGTGTCAGGGTGACCATGGTGAGCAAGGGCGAGGA", ssdna=True, product="RS204") #Create a QUEEN object for the forward primer. 
SI627 = QUEEN(seq="CTCCCCTACCCGGTAGAATTGGATCCTTACTTGTACAGCTCGTCC", ssdna=True, product="SI627")        #Create a QUEEN object for the reverse primer. 
pLV_eGFP.searchsequence(query=RS204.seq[-18:], product="FW", pn=pn, pd=pd)         #Search for the 18-bp 3’-end sequences of the forward primer.
pLV_eGFP.searchsequence(query=SI627.seq[-18:], product="RV", pn=pn, pd=pd)         #Search for the 18-bp 3’-end sequences of the reverse primer.
extract1    = cropdna(pLV_eGFP, FW[0].end, RV[0].start, product="extract1", pn=pn, pd=pd)        #Crop the internal DNA sequence flanked by the primer annealing sites.
PCR_product = modifyends(extract1, RS204.seq, SI627.rcseq, product="PCR_product", pn=pn, pd=pd)  #Add forward and reverse primer sequences to the both ends of the cropped fragment. 

pn = "Restriction enzyme digestion" 
pd = "2. The eGFP-encoding cassette was digested by EcoRI and BamHI."
EcoRI_site_ins = PCR_product.searchsequence(cs.lib["EcoRI"], product="EcoRI_site_ins", pn=pn, pd=pd)   #Search for EcoRI site.
BamHI_site_ins = PCR_product.searchsequence(cs.lib["BamHI"], product="BamHI_site_ins", pn=pn, pd=pd)   #Search for BamHI site.
fragment_1     = cropdna(PCR_product, EcoRI_site_ins[0], BamHI_site_ins[0], product="fragment_1", pn=pn, pd=pd) #Cut "PCR_product" at the cut sites.

pn = "Restriction enzyme digestion" 
pd = "3. The destination plasmid pLVSIN-CMV-Puro (Takara) was also digested by EcoRI and BamHI."
EcoRI_site_bk = pLVSIN_CMV_pur.searchsequence(cs.lib["EcoRI"], product="EcoRI_bk", pn=pn, pd=pd)   #Search for EcoRI site.
BamHI_site_bk = LVSIN_CMV_pur.searchsequence(cs.lib["BamHI"],  product="BamHI_bk", pn=pn, pd=pd)   #Search for EcoRI site.
fragment_2    = cropdna(pLVSIN_CMV_pur, BamHI_site_bk[0], EcoRI_site_bk[0], product="fragment_2", pn=pn, pd=pd) #Cut pLVSIN_CMV_pur at the cut sites.

pn     = "Ligation"
pd     = "4. The eGFP-encoding cassette was cloned into the destination vector using T4 DNA ligase (NEB)."
pRS112 = joindna(fragment_1,fragment_2, product="pRS112", topology="circular", pn=pn, pd=pd) #Join "fragment_1" and "framgnet_2" and generate "pRS112" plasmid object.
```



## Example 2: pCMV-Target-ACE  construction  

QUEEN script for pCMV-Target-ACE construction by Gibson Assembly process.

1. A backbone fragment was amplified from pCMV-ABE7.10 using the primer pair RS047/RS052.
2. The C-terminus region of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) using the primer pair RS051/RS046.
3. The two fragments were assembled by Gibson Assembly.

```python
#Load pCMV-ABE plasmid object from a Benchling share link.
pCMV_ABE = QUEEN(record="https://benchling.com/s/seq-K4HkSd2E8WiTAulJUeBf", dbtype="benchling", product="pCMV_ABE")            

#Load pcDNA31_Target_AID plasmid object from a Benchling share link.
pcDNA31_Target_AID = QUEEN(record="https://benchling.com/s/seq-cfnGDU0Mq8cUwn185LPF", dbtype="benchling", product="pcDNA31_Target_AID")  

processname1 = "PCR"
description1 = "1. A backbone fragment was amplified from pCMV-ABE7.10 using the primer set RS047/RS052."
RS052 = QUEEN(seq="ACCTCCTCCACCGTCACCCCCAAGCTGTGACA", product="RS052") #Create a QUEEN object for the forward primer. 
RS047 = QUEEN("ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC",  product="RS047") #Create a QUEEN object for the reverse primer. 
FW1   = pCMV_ABE.searchsequence(RS047.seq[-18:], product="FW1", pn=processname1, pd=description1) #Search for the 18-bp 3’-end sequences of the forward primer.
RV1   = pCMV_ABE.searchsequence(RS052.seq[-18:], product="RV1", pn=processname1, pd=description1) #Search for the 18-bp 3’-end sequences of the reverse primer.
extract1  = cropdna(pCMV_ABE, FW4[0].end, RV4[0].start, product="extract1", pn=processname1, pd=description1)   #Crop the internal DNA sequence flanked by the primer annealing sites.
fragment1 = modifyends(extract4, RS047.seq, RS052.rcseq, product="fragment1", pn=processname1, pd=description1) #Add forward and reverse primer sequences to the both ends of the cropped fragment. 

processname2 = "PCR"
description2 = "2. The C-terminus region of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using the primer set RS051/RS046."
RS046 = QUEEN(seq="TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC", product="RS046") #Create a QUEEN object for the forward primer. 
RS051 = QUEEN(seq="GCTTGGGGGTGACGGTGGAGGAGGTACCGGCGG",    product="RS051") #Create a QUEEN object for the reverse primer. 
FW2   = pcDNA31_Target_AID.searchsequence(RS051.seq[-18:], product="FW2", pn=processname2, pd=description2) #Search for the 18-bp 3’-end sequences of the forward primer.
RV2   = pcDNA31_Target_AID.searchsequence(RS046.seq[-18:], product="RV2", pn=processname2, pd=description2) #Search for the 18-bp 3’-end sequences of the reverse primer.
extract2  = cropdna(pcDNA31_Target_AID, FW5[0].end, RV5[0].start, product="extract2", pn=processname2, pd=description2) #Crop the internal DNA sequence flanked by the primer annealing sites.
fragment2 = modifyends(extract5, RS051.seq, RS046.rcseq, product="fragment2", pn=processname2, pd=description2)         #Add forward and reverse primer sequences to the both ends of the cropped fragment. 

processname3 = "Gibson Assembly"
description3 = "3. The two fragments were assembled by Gibson Assembly reaction."
fragment1 = modifyends(fragment1, "*{25}/-{25}","-{25}/*{25}", product="fragment1", pn=processname3, pd=description3) #Generate long sticky ends on the both sides of "fragment1".
fragment2 = modifyends(fragment2, "*{25}/-{25}","-{25}/*{25}", product="fragment2", pn=processname3, pd=description3) #Generate long sticky ends on the both sides of "fragment2".
pCMV_Target_ACE = joindna(fragment1, fragment2, topology="circular", product="pCMV_Target_ACE", pn=processname3, pd=description3) #Join the fragments.
```

