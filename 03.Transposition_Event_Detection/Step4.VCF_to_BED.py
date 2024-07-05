from Bio import Seq
from Bio import SeqIO
import os
import pandas as pd
def get_unique_sv(sv_df):
    sv_df=sv_df.reset_index(drop=True)
    pos_samples=[]
    for s in list(sv_df.columns)[9:-5]:
        temp=sv_df[(sv_df[s]!='./.')&(sv_df[s]!='0/0')]
        if temp.shape[0]>0:
            pos_samples.append(s)
    sv_df=sv_df.drop_duplicates("SV")
    sv_df.loc[0,pos_samples]="0/1"
    sv_df["Support_SVs"]=len(pos_samples)
    return sv_df
if 'SV_Bed' not in os.listdir("/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/VCF/Brassica_oleracea/"):
    os.mkdir("/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/VCF/Brassica_oleracea/SV_Bed/")
if 'Merged_VCF_Less' not in os.listdir("/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/VCF/Brassica_oleracea/"):
    os.mkdir("/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/VCF/Brassica_oleracea/Merged_VCF_Less/")
for file in os.listdir('/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/VCF/Brassica_oleracea/Merged_VCF'):
    if ".vcf" in file:
        ref_sample=file.split('.1000bp.merged')[0]
        if file.endswith(".gz"):
            df=pd.read_csv('/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/VCF/Brassica_oleracea/Merged_VCF/'+file,compression='gzip',sep='\t',skiprows=33)
        elif file.endswith(".vcf"):
            df=pd.read_csv('/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/VCF/Brassica_oleracea/Merged_VCF/'+file,sep='\t',skiprows=33)
        columns=df.columns
        samples=list(df.columns[9:])
        for sample in samples:
            df[sample]=df[sample].apply(lambda x:x.split(":")[0])
        df['FORMAT']='GT'
        df["SV_Type"]=df['INFO'].apply(lambda x:x.split("SVTYPE=")[1].split(";")[0])   
        df["SV"]=df["#CHROM"].apply(lambda x:str(x))+":"+df["POS"].apply(lambda x:str(x))+":"+df["SV_Type"]
        df["Support_SVs"]=df['INFO'].apply(lambda x:int(x.split("SUPP=")[1].split(";")[0]))
        df["SV_Length"]=df['INFO'].apply(lambda x:abs(int(x.split("SVLEN=")[1].split(";")[0])))
        df["nonN_Percent"]=(df["REF"].apply(lambda x:len(str(x))-str(x).upper().count("N"))+df["ALT"].apply(lambda x:len(str(x))-str(x).upper().count("N")))/(df["REF"].apply(lambda x:len(str(x)))+df["ALT"].apply(lambda x:len(str(x))))
        df["INFO"]=df["INFO"].apply(lambda x:x.split(";SVMETHOD")[0])
		
        new=[]
        count=0
        for sv,sv_df in df.groupby("SV"):
            sv_df=sv_df.sort_values("nonN_Percent",ascending=False).reset_index(drop=True)
            count+=1
            if count%5000==0:
                print(count)
            if sv_df.shape[0]==1:
                new.append(sv_df)
            elif sv_df.shape[0]>1:
                sv_df=get_unique_sv(sv_df)
                new.append(sv_df)
        new=pd.concat(new)
	#print(new.shape)
        df=new
        df["ID"]=df["SV"]
        df["End"]=df["POS"]+df["REF"].apply(lambda x:len(str(x)))-1
		
        bed=df[['#CHROM',"POS","End","ID",'Support_SVs']]
        bed.columns=['#CHROM',"Start","End","ID",'Support_SVs']
        bed['#CHROM']=bed["#CHROM"].apply(lambda x:x.replace(ref_sample+".",''))
        bed.to_csv("/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/VCF/Brassica_oleracea/SV_Bed/"+ref_sample+"_SV.bed",sep='\t',index=False)
		
        df["INFO"]='SUPP='+df['Support_SVs'].apply(lambda x:str(x))+";SUPP_VEC="+df["INFO"].apply(lambda x:x.split(";SUPP_VEC=")[1])
        vcf=df[columns]
        vcf=vcf.sort_values(["#CHROM","POS"])
        vcf.to_csv('/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/VCF/Brassica_oleracea/Merged_VCF_Less/'+ref_sample+'.sv.vcf',sep='\t',index=False)
		
