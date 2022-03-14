import numpy as np
import pandas as pd
import pymysql
from math import log
import math


def OR_Calculate(user_allele, risk_allele):
    if risk_allele == '?':
        return 0
    else:
        add_score = 0
        add_score = user_allele.count(risk_allele)
        return add_score

def MAF_Calculate(user_allele, risk_allele):
    if risk_allele == '?':
        return 0
    else:
        MAF=0
        MAF=user_allele
        weighted_MAF=math.sqrt(2*MAF*(1-MAF))
        return weighted_MAF


def Method3_MAF(Disease_name):
    ############################   Step 1    #########################################
    print("{:=^50s}".format("1.数据准备、匹配SNP"))
    # 1.1 导入用户基因数据
    # 处理数据，删除无效信息，按照rs_ID排序(方便后续匹配查找)
    df_user = pd.read_csv('snp-inf.csv',
                          usecols=['dbSNP_RS_ID', '1.CEL_call_code', '2.CEL_call_code', '3.CEL_call_code', 'Minor_Allele_Frequency'])
    user=df_user.loc[:,['dbSNP_RS_ID','Minor_Allele_Frequency']]
    user=user.reset_index(drop=True)
    valid_df_user = df_user[df_user['dbSNP_RS_ID'] != 'nan']
    sorted_valid_df_user = valid_df_user.sort_values('dbSNP_RS_ID')
    user_rs = sorted_valid_df_user['dbSNP_RS_ID']
    user_rs_num = len(user_rs)
    print('有效提供的SNP位点数量:  ', user_rs_num)

    # 1.2 导入GWAS
    simplified_database_GWAS = pd.read_csv('simplified_database_GWAS.csv',
                                           usecols=['English_GWAS', 'Risk_allele', 'OR'])
    Specific_GWAS = simplified_database_GWAS[simplified_database_GWAS['English_GWAS'] == Disease_name].loc[:,
                    ['Risk_allele', 'OR']]  # 特定疾病在GWAS库里有的SNP以及对应OR值
    Specific_GWAS = Specific_GWAS.reset_index(drop=True)  # 重置index
    GWAS_rs = sorted(set([str(x).split('-')[0] for x in Specific_GWAS['Risk_allele']]))  # 提取GWAS中的SNP（去重）
    print('在GWAS库中该疾病有效的风险SNP个数： ', len(Specific_GWAS))
    if len(Specific_GWAS) == 0:
        # 风险得分，GWAS中有效SNP个数，可匹配SNP个数，SNP相关信息 or 情况说明，是否无信息
        return [[0, 0, 0], 0, 0, 'GWAS数据库中无SNP信息', 0]

        # 1.3 匹配GWAS中与提供信息文件中重合的SNP
    user_rs_index = []  # 记录匹配上SNP在user表中的序号
    matched_rs = []  # 记录匹配上的SNP
    for rs in GWAS_rs:
        for index, user_rs in zip(range(len(df_user['dbSNP_RS_ID'])), df_user['dbSNP_RS_ID']):
            if rs == user_rs:
                matched_rs.append(rs)
                user_rs_index.append(index)

    Num_Match_rs = len(matched_rs)
    print("能匹配上的SNP个数：  ", Num_Match_rs)
    if Num_Match_rs == 0:
        # 风险得分，GWAS中有效SNP个数，可匹配SNP个数，SNP相关信息 or 情况说明，是否无信息
        return [[0, 0, 0], len(Specific_GWAS), Num_Match_rs, '无匹配上的SNP', 0]
    print(matched_rs)  # 能配对上的SNP

    ############################   Step 2    #########################################
    print("{:=^50s}".format("2.Calculate Weights"))
    # 2.1 建立字典Dict存储每个SNP各风险等位基因的信息
    Dict = dict()
    # 风险基因&OR值存入字典Dict(不包含用户不提供的SNP位点信息)
    # Dict={'rs9592461': [['A', [7.062, 1.025]]], ['G', [1.021451]]], 'rs645592': [['T', [5.571]]]}
    for index, row in Specific_GWAS.iterrows():
        rs = str(row['Risk_allele']).split('-')[0]
        AGCT = str(row['Risk_allele']).split('-')[1][-1]
        if rs in matched_rs:
            if rs in Dict:  # rs在字典中存在
                if AGCT not in [i[0] for i in Dict[rs]]:  # 但AGCT尚未存在时
                    Dict[rs].append([AGCT, [float(row['OR'])]])
                else:
                    for i in Dict[rs]:
                        if i[0] == AGCT:
                            i[1].append(float(row['OR']))
            elif not (row['OR'] is None):
                Dict[rs] = [[AGCT, [float(row['OR'])]]]
    # print(Dict)


    # 2.2 根据字典Dict计算各SNP等位风险基因对应的风险权重
    for key, value in Dict.items():
        # value,eg.:[['A', [7.062, 1.025]]], ['G', [1.021451]]]
        n = len(value)
        if n > 0:
            for x in value:
                # x,eg.:['A', [7.062, 1.025]]
                AGCT = x[0]
                OR_list = x[1]
                OR_num = len(x[1])
                OR = sum(OR_list) / OR_num
                OR_log = log(OR)
                x[1] = OR_log
                # x,eg.:['A', 0.3]
            Dict[key] = value
    print('使用到的位点信息： ', Dict)  # {'rs1656369': [['A', 1.76541],['G', 1.76541549]], 'rs10811965': [['T', 1.7279318]]}

    ##############################       Step 3    ###########################################
    print("{:=^50s}".format("3.Calculate Results"))
    # 3.1 计算risk_max
    risk_max = 0
    for rs in Dict:
        for i in Dict[rs]:
            if i[1] != 0:
                risk_max += i[1] * 2
    risk_max = abs(risk_max)  # 取绝对值
    print('risk_max   ', risk_max)

    # 3.2 计算分数
    risk_score_1, risk_score_2, risk_score_3 = 0, 0, 0
    for index in user_rs_index:
        rs = sorted_valid_df_user['dbSNP_RS_ID'][index]
        if rs in Dict:
            for i in Dict[rs]:
                risk_allele = i[0]
                weighted_OR = i[1]
                #weighted_MAF=math.sqrt(i[3]*2*(1-i[3]))
                MAF=sorted_valid_df_user['Minor_Allele_Frequency'][index]
                MAF1=MAF[0:7]
                MAF2=MAF[18:25]
                MAF3=MAF[36:43]
                MAF4=MAF[54:61]
                MAF0=float((float(MAF1)+float(MAF2)+float(MAF3)+float(MAF4))/4)
                weighted_MAF=math.sqrt(2*float(MAF0)*(1-float(MAF0)))
                risk_score_1 += weighted_MAF*weighted_OR * OR_Calculate(sorted_valid_df_user['1.CEL_call_code'][index], risk_allele)
                risk_score_2 += weighted_MAF*weighted_OR * OR_Calculate(sorted_valid_df_user['2.CEL_call_code'][index], risk_allele)
                risk_score_3 += weighted_MAF*weighted_OR * OR_Calculate(sorted_valid_df_user['3.CEL_call_code'][index], risk_allele)

    # 归一化结果
    final_score = [risk_score_1 / risk_max, risk_score_2 / risk_max, risk_score_3 / risk_max]
    print("三位的患病风险： ", final_score)
    # 风险得分，GWAS中有效SNP个数，可匹配SNP个数，SNP相关信息 or 情况说明，是否无信息
    return [final_score, len(Specific_GWAS), Num_Match_rs, Dict, 1]


