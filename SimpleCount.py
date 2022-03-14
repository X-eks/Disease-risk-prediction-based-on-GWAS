import numpy as np
import pandas as pd
import vcf2excel

def OR_Calculate(user_allele,risk_allele):
    if risk_allele =='?':
        return 0
    else:
        add_score=0
        add_score=user_allele.count(risk_allele)
        return add_score




#简单相加遗传风险评分（a simple count genetic risk score，SC-GRS）
def Method2_SC(Dict, df_user, user_rs_index):
    
    #计算risk_max
    risk_max = 0
    for rs in Dict:
        risk_max += len(Dict[rs])
    risk_max=risk_max*2
    #print('risk_max   ',risk_max)

    #计算分数
    risk_score=0
    for index in user_rs_index:
        rs=df_user['rs_ID'][index]
        if rs in Dict:
            for i in Dict[rs]:
                risk_allele = i[0]        
                risk_score += OR_Calculate(df_user['Allele_1'][index],risk_allele)
                risk_score += OR_Calculate(df_user['Allele_2'][index],risk_allele)

    return risk_score/risk_max

