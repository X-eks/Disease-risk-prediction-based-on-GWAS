import pandas as pd
    
    
#风险对应表--1（文字版）
def Risk_Level_1(RR):
    if RR<0.5:
        return '显著低于一般人群'
    elif RR<0.9:
        return '低于一般人群'
    elif RR<1.1:
        return '相当于一般人群'
    elif RR<2.0:
        return '高于一般人群'
    else:
        return '显著高于一般人群'

def score2level(Disease_name, score, typee):  

    RR_avg=pd.read_excel('RR_avg.xlsx')
    for index, row in RR_avg.iterrows():
        if row['English']== Disease_name:

            if typee=='OR':
                if row['OR_count']==0:
                    return Risk_Level_1(1)
                elif row['OR_score_sum']==0:
                    return '群体风险值为0'
                else:
                    group_avg= row['OR_score_sum']/row['OR_count']

            elif typee=='SC':
                if row['SC_count']==0:
                    return Risk_Level_1(1)
                elif row['SC_score_sum']==0:
                    return '群体风险值为0'
                else:
                    group_avg= row['SC_score_sum']/row['SC_count']

            elif typee=='SD':
                if row['SD_count']==0:
                    return Risk_Level_1(1)
                elif row['SD_score_sum']==0:
                    return '群体风险值为0'
                else:
                    group_avg= row['SD_score_sum']/row['SD_count']
            break

    print('Group_avg:  ', group_avg)
    risk_level=Risk_Level_1(score/group_avg)
    return risk_level