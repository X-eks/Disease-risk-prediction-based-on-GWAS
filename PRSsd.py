import numpy as np
import pandas as pd






def Method3_PRSsd(disease_GWAS, OR_score):
    # 更新表格数据
    risk_score=0
    flag=False# 是否插入标志
    df_SD_score=pd.read_excel('SD_score.xlsx')
    for index, row in df_SD_score.iterrows():
        if row[0]==disease_GWAS:
            df1=df_SD_score.iloc[:index,:]
            df2=df_SD_score.iloc[index:,:]
            df_new=pd.concat([df1, pd.Series([disease_GWAS, OR_score]),df2], ignore_index=True)
            break
    if flag==False:
        df_new=df_SD_score.append({'English': disease_GWAS, 'OR_score': OR_score}, ignore_index=True)
    df_new.to_excel('SD_score.xlsx', index=False)
    #print(df_new)


    #计算均值
    li=[]
    for index, row in df_new.iterrows():
        if row[0]==disease_GWAS:
            li.append(float(row[1]))
    print('PRSsd表格中的数据：',li)

    prsmean=np.mean(li)
    prssd=np.std(li)
    if len(li)==0 or prssd==0:
        risk_score=OR_score
    else:
        risk_score=abs(risk_score-prsmean)/prssd

    return risk_score





