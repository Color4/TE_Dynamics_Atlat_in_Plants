from Bio import Seq
from Bio import SeqIO
import os
import pandas as pd
def extract_Ref_sequences(loc):
    chrom=loc.split(":")[0]
    start = int(loc.split(":")[1].split("-")[0]) - 1  # 转换为 0-based 索引
    end = int(loc.split(":")[1].split("-")[1])
    strand = loc.split(":")[2]
    if chrom in ref_dic:
        if strand == "-":
            seq = ref_dic[chrom][end:start]   
            seq = Seq.reverse_complement(seq) # 考虑到负链的情况，需要进行反向互补       
        else:
            seq = ref_dic[chrom][start:end] 
        return seq
    else:
        seq=""
        return seq
def extract_Alt_sequences(loc):
    chrom=loc.split(":")[0]
    start = int(loc.split(":")[1].split("-")[0]) - 1  # 转换为 0-based 索引
    end = int(loc.split(":")[1].split("-")[1])
    strand = loc.split(":")[2]
    if chrom in query_dic:
        if strand == "-":
            seq = query_dic[chrom][end:start]   
            seq = Seq.reverse_complement(seq) # 考虑到负链的情况，需要进行反向互补       
        else:
            seq = query_dic[chrom][start:end] 
        return seq
    else:
        seq=""
        return seq
path='/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/'
species=["Brassica_oleracea"]#,,'Maize','Tomato',"Maize"
for specie in species:
    if 'Raw' not in os.listdir(path+"/VCF/"+specie):
        os.mkdir(path+"/VCF/"+specie+"/Raw/")
    file_path=path+specie+"/"
    R=[]
    samples=[sample.split("_SV")[0] for sample in os.listdir(file_path)]
    for ref in samples:
        print(ref)
        #ref_dic=fa_dic[ref]
        ref_dic={}
        for s in SeqIO.parse("/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/TE_Mobile/Genome/"+specie+"/"+ref+'.fa','fasta'):
            ref_dic[s.id]=str(s.seq)
        results=os.listdir(file_path+ref+"_SV")
        if ref not in os.listdir(path+"/VCF/"+specie+"/Raw/"):
            os.mkdir(path+"/VCF/"+specie+"/Raw/"+ref)
        for result in results:
            #print(df.shape)
            lines= open(file_path+ref+"_SV/"+result,'r').readlines()
            if len(lines)>5:
                df=pd.read_csv(file_path+ref+"_SV/"+result,sep='\t').dropna()
                df=df[(df["LEN"].apply(lambda x:abs(x))<1000000)&(df["LEN"].apply(lambda x:abs(x))>=50)].reset_index(drop=True)
                df["Ref_Strand"]="+"
                idxs=df[df["REF_START"]>df["REF_END"]].index
                df.loc[idxs,"Ref_Strand"]="-"
                df["Alt_Strand"]="+"
                idxs=df[df["Q_START"]>df["Q_END"]].index
                df.loc[idxs,"Alt_Strand"]="-"
                query=result.split("sv.")[1].split("_"+ref)[0]
                #if ref in ["Zm-B73",'Zm-SK']:
                #    df["Q_CHROM"]=df["Q_CHROM"].apply(lambda x:x.replace(query+".",""))
                    #df["Alt_Loc"]=df["Q_CHROM"].apply(lambda x:x.replace(query+".",""))+":"+df["Q_START"].apply(lambda x:str(x))+"-"+df["Q_END"].apply(lambda x:str(x))+":"+df["Alt_Strand"]
                    #df["Ref_Loc"]=df["REF_CHROM"]+":"+df["REF_START"].apply(lambda x:str(x))+"-"+df["REF_END"].apply(lambda x:str(x))+":"+df["Ref_Strand"]
                #else:
                #    df["REF_CHROM"]=df["REF_CHROM"].apply(lambda x:x.replace(ref+".",""))
                    #df["Ref_Loc"]=df["REF_CHROM"].apply(lambda x:x.replace(query+".",""))+":"+df["REF_START"].apply(lambda x:str(x))+"-"+df["REF_END"].apply(lambda x:str(x))+":"+df["Ref_Strand"]
                df["Ref_Loc"]=df["REF_CHROM"]+":"+df["REF_START"].apply(lambda x:str(x))+"-"+df["REF_END"].apply(lambda x:str(x))+":"+df["Ref_Strand"]
                df["Alt_Loc"]=df["Q_CHROM"]+":"+df["Q_START"].apply(lambda x:str(x))+"-"+df["Q_END"].apply(lambda x:str(x))+":"+df["Alt_Strand"]
                query_dic={}
                for s in SeqIO.parse("/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/TE_Mobile/Genome/"+specie+"/"+query+'.fa','fasta'):
                    query_dic[s.id]=str(s.seq)
                df["Ref_Seq"]=df["Ref_Loc"].apply(extract_Ref_sequences)
                df["Alt_Seq"]=df["Alt_Loc"].apply(extract_Alt_sequences)
                df=df[(df["Ref_Seq"]!="")&(df["Alt_Seq"]!="")].reset_index(drop=True)
                df=df.sort_values(["REF_CHROM","REF_START"])
                df["Alt_N_Ratio"]=df["Alt_Seq"].apply(lambda x:round(str(x).count("N")/len(str(x)),2))
                df["Ref_N_Ratio"]=df["Ref_Seq"].apply(lambda x:round(str(x).count("N")/len(str(x)),2))
                df=df[((df["SV_TYPE"]=="INS")&(df["Alt_N_Ratio"]<0.9))|
                  ((df["SV_TYPE"]!="INS")&(df["SV_TYPE"]!="DEL"))|
                  ((df["SV_TYPE"]=="DEL")&(df["Ref_N_Ratio"]<0.9))]
                df["ID"]=df["Ref_Loc"]
                df["XI"]=df["Alt_Loc"]
                df["INFO"]="SVTYPE="+df["SV_TYPE"]+";NA=1;XI="+df["XI"]+";CE="+df["SV_TYPE"]
                df["QUAL"]="."
                df["FILTER"]="."
                df["FORMAT"]="GT"
                df[query]='0/1'
                df=df[['REF_CHROM', 'REF_START', 'ID', 'Ref_Seq', 'Alt_Seq', "QUAL","FILTER", 'INFO',"FORMAT",query]]
                df.columns=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",query]
                print(query,df.shape)
                df.to_csv(path+"/VCF/"+specie+"/Raw/"+ref+"/"+"sv."+query+"_"+ref+".vcf",sep='\t',index=False)
                del df
