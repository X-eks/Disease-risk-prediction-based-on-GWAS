import pandas as pd
import re
import os

def vcf_process(file_dir):
    with open(file_dir, 'r', encoding='utf-8') as rf:
        column = []
        datas = []
        #提取file_ID
        extract_temp= re.findall(r'(\d+|\D+)', file_dir)
        numbers=[]
        for item in extract_temp:
            if item.isdigit():
                numbers.append(item)
        file_ID=numbers[-1]

        for line in rf.readlines():
            if line.startswith('#') and line.startswith("##") == False: #跳过##部分
                column=[ i for i in str(line).lstrip('#').split()] #列标签
                #print(column)
            elif line.startswith('#') == False:
                datas.append([ i for i in str(line).split()]) #数据      
        data = pd.DataFrame(datas,columns = column)#从.vcf中取出的数据

    #创建新文件
    a = []
    b = []
    c = []
    for index, row in data.iterrows():
        if index%100000==0:
            print('格式转换进度：  ', str(float(index)/float(data.shape[0])*100)[:5], '%')
        try:
            rs_info=re.findall(r'RSID=[\w]*',row['INFO'])[0]   #截取rsID
            rs_ID=rs_info[5:]
            
            abbr_code=row['s'+ file_ID]
            if abbr_code[0] == '0' and abbr_code[2] == '0': #根据s1得到REF和ALT
                b.append(row['REF']) 
                c.append(row['REF'])
                a.append(rs_ID)
            elif abbr_code[0] == '0' and abbr_code[2] == '1':
                b.append(row['REF'])
                c.append(row['ALT'])
                a.append(rs_ID)
            elif abbr_code[0] == '1' and abbr_code[2] == '1':
                b.append(row['ALT'])
                c.append(row['ALT'])
                a.append(rs_ID)
        except Exception as e:
            pass
            #print(e)   #输出错误信息的具体情况
    data_dict={'rs_ID': a, 'Allele_1': b,'Allele_2': c} 
    user_data = pd.DataFrame(data_dict)
    return user_data