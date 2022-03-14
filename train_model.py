import SimpleCount
import OddsRatio
import PRSsd
import pandas as pd
import numpy as np
import vcf2excel
from math import log
import score2level




def train_model(file_dir):   #输入被测试样本文件的路径

    #创建输出表格
    Valid_Results=pd.DataFrame(columns=['Chinese', 'English', 'OR_score','SC_score','SD_score', 'SNP_info'])
    Report_data=pd.DataFrame(columns=['Chinese', 'English', 'OR_level','SC_level','SD_level', 'SNP_info'])
    SNP_info=[]
    Num_Match_SNP=[]

    #导入用户基因数据
    df_user=vcf2excel.vcf_process(file_dir)
    user_rs=df_user['rs_ID']#按照rs_ID排序(方便后续匹配查找)
    user_rs_num = len(user_rs)
    print('用户表里的信息总数：', user_rs_num)


    #199种带OR值的疾病数据
    ##导入数据
    Disease_Avai=pd.read_csv('199_disease_name.csv', usecols=['English_GWAS'])['English_GWAS'].tolist()
    Disease_Chinese=pd.read_csv('199_disease_name.csv', usecols=['Chinese'])['Chinese'].tolist()
    Disease_Avai_test=Disease_Avai
    ##计算分数
    print('-'*20+"  数据库1  "+'-'*20)
    for disease_GWAS, index in zip(Disease_Avai_test, range(len(Disease_Avai_test))):
        Chinese_name=Disease_Chinese[index]
        print("{:=^50s}".format('当前疾病： '+ disease_GWAS))
        print('当前运行进度： ' , index+1, '/', len(Disease_Avai_test))

        #获取GWAS数据库中关于该疾病的信息
        simplified_database_GWAS=pd.read_csv('simplified_database_GWAS.csv', usecols=['English_GWAS','Risk_allele','OR'])
        Specific_GWAS=simplified_database_GWAS[simplified_database_GWAS['English_GWAS']==disease_GWAS].loc[:, ['Risk_allele','OR']]#特定疾病在GWAS库里有的SNP以及对应OR值
        Specific_GWAS=Specific_GWAS.reset_index(drop=True)#重置index
        GWAS_rs = sorted(set([str(x).split('-')[0] for x in Specific_GWAS['Risk_allele']]))#提取GWAS中的SNP（去重）
        print('在GWAS库中该疾病有效的风险SNP个数： ', len(Specific_GWAS))
        if len(Specific_GWAS)==0:
            print('err1：GWAS数据库中无SNP信息')
            Err_info='err1'
            continue

        #匹配GWAS中与用户基因文件中重合的SNP
        user_rs_index = []  #记录匹配上SNP在user表中的序号
        matched_rs = [] #记录匹配上的SNP
        for rs in GWAS_rs:
            for index, user_rs in zip(range(len(df_user['rs_ID'])), df_user['rs_ID']):
                if rs==user_rs:
                    matched_rs.append(rs)
                    user_rs_index.append(index)

        Num_Match_rs=len(matched_rs)
        print("能匹配上的SNP个数：  ", Num_Match_rs)
        if Num_Match_rs==0:
            print('err2：无匹配上的SNP')
            Err_info='err2'
            continue
        print(matched_rs)  # 能配对上的SNP

        #建立权重字典：Dict_OR & Dict_SC （建立字典Dict存储每个SNP各风险等位基因的信息）(目前三种模型可以使用)
        #####   Dict_OR
        Dict_OR = dict() 
        #风险基因&OR值存入字典Dict(不包含用户不提供的SNP位点信息)
        #Dict_OR={'rs9592461': [['A', [7.062, 1.025]]], ['G', [1.021451]]], 'rs645592': [['T', [5.571]]]}
        for index, row in Specific_GWAS.iterrows():
            rs = str(row['Risk_allele']).split('-')[0]
            AGCT = str(row['Risk_allele']).split('-')[1][-1]
            if rs in matched_rs:
                if rs in Dict_OR : #rs在字典中存在
                    if AGCT not in [i[0] for i in Dict_OR[rs]]:  #但AGCT尚未存在时
                        Dict_OR[rs].append([AGCT,[float(row['OR'])]])
                    else:
                        for i in Dict_OR[rs]:
                            if i[0]==AGCT:
                                i[1].append(float(row['OR']))
                elif not(row['OR'] is None):
                    Dict_OR[rs]=[[AGCT,[float(row['OR'])]]]

        #根据字典Dict计算各SNP等位风险基因对应的风险权重
        for key, value in Dict_OR.items():
            #value,eg.:[['A', [7.062, 1.025]]], ['G', [1.021451]]]
            n = len(value)
            if n > 0:
                for x in value:
                    #x,eg.:['A', [7.062, 1.025]]
                    AGCT = x[0]
                    OR_list = x[1]
                    OR_num=len(x[1])
                    OR=sum(OR_list)/OR_num
                    OR_log=log(OR)
                    x[1]=OR_log
                    #x,eg.:['A', 0.3]
                Dict_OR[key]=value
        print('使用到的位点信息： ', Dict_OR)  #{'rs1656369': [['A', 1.76541],['G', 1.76541549]], 'rs10811965': [['T', 1.7279318]]}


        #####   Dict_SC
        Dict_SC = dict() 
        for index, row in Specific_GWAS.iterrows():
            rs = str(row['Risk_allele']).split('-')[0]
            AGCT = str(row['Risk_allele']).split('-')[1][-1]
            if rs in matched_rs and not(row['OR'] is None):
                if rs in Dict_SC and AGCT not in Dict_SC[rs]:
                    Dict_SC[rs].append(AGCT)
                else:
                    Dict_SC[rs]=[AGCT] 
        print(Dict_SC) #{'rs1234:['A','G'], 'rs1235':['C']}


        #运行模型
        OR_score= OddsRatio.Method1_OR(Dict_OR, df_user, user_rs_index)
        SC_score= SimpleCount.Method2_SC(Dict_SC, df_user, user_rs_index)
        SD_score= PRSsd.Method3_PRSsd(disease_GWAS, OR_score)

        #更新分数记录表格
        update_RR_avg(True, [OR_score, SC_score, SD_score], disease_GWAS)
        Valid_Results=Valid_Results.append({'Chinese': Chinese_name, 'English': disease_GWAS,
        'OR_score': OR_score, 'SC_score': SC_score, 'SD_score': SD_score, 'SNP_info': Dict_OR}, ignore_index=True)
        Report_data=Report_data.append({'Chinese': Chinese_name, 'English': disease_GWAS,
        'OR_level': score2level.score2level(disease_GWAS, OR_score, 'OR'), 'SC_level': score2level.score2level(disease_GWAS, SC_score, 'SC'),\
             'SD_level': score2level.score2level(disease_GWAS, SD_score, 'SD'), 'SNP_info': Dict_OR}, ignore_index=True)




    #导入35种不带OR值的疾病数据（目前只有SC一种模型）
    print('-'*20+"  数据库2  "+'-'*20)
    ##导入数据
    Disease_Avai=pd.read_excel('35_disease_name.xlsx')['English'].tolist()
    Disease_Chinese=pd.read_excel('35_disease_name.xlsx')['Chinese'].tolist()
    Disease_Avai_test=Disease_Avai

    for disease_GWAS, index in zip(Disease_Avai_test, range(len(Disease_Avai_test))):
        print('当前疾病： ', disease_GWAS)
        print('当前运行进度： ' , index+1, '/', len(Disease_Avai_test))
        Chinese_name=Disease_Chinese[index]

        #获取GWAS数据库中关于该疾病的信息
        No_OR_database=pd.read_excel('35_disease_database.xlsx')
        No_OR_database=No_OR_database[No_OR_database['English']==disease_GWAS].loc[:, ['Risk_allele']]#特定疾病在GWAS库里有的SNP以及对应OR值
        No_OR_database=No_OR_database.reset_index(drop=True)#重置index
        GWAS_rs = sorted(set([str(x).split('-')[0] for x in No_OR_database['Risk_allele']]))#提取GWAS中的SNP（去重）
        print('在GWAS库中该疾病有效的风险SNP个数： ', len(No_OR_database))
        if len(No_OR_database)==0:
            print('err1：GWAS数据库中无SNP信息')
            Err_info='err1'
            continue

        #匹配GWAS中与用户基因文件中重合的SNP
        user_rs_index = []  #记录匹配上SNP在user表中的序号
        matched_rs = [] #记录匹配上的SNP
        for rs in GWAS_rs:
            for index, user_rs in zip(range(len(df_user['rs_ID'])), df_user['rs_ID']):
                if rs==user_rs:
                    matched_rs.append(rs)
                    user_rs_index.append(index)

        Num_Match_rs=len(matched_rs)
        print("能匹配上的SNP个数：  ", Num_Match_rs)
        if Num_Match_rs==0:
            print('err2：无匹配上的SNP')
            Err_info='err2'
            continue
        print(matched_rs)  # 能配对上的SNP
    
        #####   Dict_SC
        Dict_SC = dict() 
        for index, row in No_OR_database.iterrows():
            rs = str(row['Risk_allele']).split('-')[0]
            AGCT = str(row['Risk_allele']).split('-')[1][-1]
            if rs in matched_rs :
                if rs in Dict_SC and AGCT not in Dict_SC[rs]:
                    Dict_SC[rs].append(AGCT)
                else:
                    Dict_SC[rs]=[AGCT] 
        print(Dict_SC) #{'rs1234:['A','G'], 'rs1235':['C']}

        #计算风险分数
        SC_score= SimpleCount.Method2_SC(Dict_SC, df_user, user_rs_index)
        #更新分数记录表格
        update_RR_avg(False, SC_score, disease_GWAS)       
        Valid_Results=Valid_Results.append({'Chinese': Disease_Chinese[index], 'English': disease_GWAS,
        'SC_score': SC_score, 'SNP_info': Dict_SC} , ignore_index=True)
        Report_data=Report_data.append({'Chinese': Chinese_name, 'English': disease_GWAS,
        'SC_level': score2level.score2level(disease_GWAS, SC_score, 'SC'), 'SNP_info': Dict_OR}, ignore_index=True)


    #将结果整理并输出至文件
    file_prefix=file_dir.split("/")[-1]
    Valid_Results.to_excel('output_data/valid/'+file_prefix+'_valid.xlsx',index=False)
    Report_data.to_excel('output_data/report_data/'+file_prefix+'_report_data.xlsx',index=False)

    

def update_RR_avg(typee, update_info, disease_GWAS):
    RR_avg=pd.read_excel('RR_avg.xlsx')
    print('RR_avg中需要更新的数据：', update_info)

    if typee==True: #199种疾病的情况，更新全部列
        for index, row in RR_avg.iterrows():
            if row[0]==disease_GWAS:
                if np.isnan(row['OR_score_sum']):
                    RR_avg.loc[index, 'OR_score_sum']=float(update_info[0])
                else:
                    RR_avg.loc[index, 'OR_score_sum']=RR_avg.loc[index, 'OR_score_sum']+float(update_info[0])
            
                if np.isnan(row['OR_count']):
                    RR_avg.loc[index, 'OR_count']=float(1)
                else:
                    RR_avg.loc[index, 'OR_count']=RR_avg.loc[index, 'OR_count']+float(1)

                if np.isnan(row['SC_score_sum']):
                    RR_avg.loc[index, 'SC_score_sum']=float(update_info[1])
                else:
                    RR_avg.loc[index, 'SC_score_sum']=RR_avg.loc[index, 'SC_score_sum']+float(update_info[1])
            
                if np.isnan(row['SC_count']):
                    RR_avg.loc[index, 'SC_count']=float(1)
                else:
                    RR_avg.loc[index, 'SC_count']=RR_avg.loc[index, 'SC_count']+float(1)

                if np.isnan(row['SD_score_sum']):
                    RR_avg.loc[index, 'SD_score_sum']=float(update_info[2])
                else:
                    RR_avg.loc[index, 'SD_score_sum']=RR_avg.loc[index, 'SD_score_sum']+float(update_info[2])
            
                if np.isnan(row['SD_count']):
                    RR_avg.loc[index, 'SD_count']=float(1)
                else:
                    RR_avg.loc[index, 'SD_count']=RR_avg.loc[index, 'SD_count']+float(1)
                print('RR_avg文件更新的这一行信息为: ', RR_avg.loc[index])
                break

    else: #35种疾病，只更新SC_score
        for index, row in RR_avg.iterrows():
            if row[0]==disease_GWAS:
                if np.isnan(row[3]):
                    RR_avg.loc[index, 'SC_score_sum']=float(update_info)
                else:
                    RR_avg.loc[index, 'SC_score_sum']=RR_avg.loc[index, 'SC_score_sum']+float(update_info)
            
                if np.isnan(row[4]):
                    RR_avg.loc[index, 'SC_count']=float(1)
                else:
                    RR_avg.loc[index, 'SC_count']=RR_avg.loc[index, 'SC_count']+float(1)
                print('RR_avg文件更新的这一行信息为: ', RR_avg.loc[index])
                break

    RR_avg.to_excel('RR_avg.xlsx')






















