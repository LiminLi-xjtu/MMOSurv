import os
import pandas as pd
import numpy as np

# 导入CSV安装包
# import csv


# cancer_list = ['BRCA', 'KIRC', 'HNSC', 'LUAD']
source_file_path = '/share/home/4120107034/MMOSurv/Datasets_and_Preprocessing/data/DNA_methylation_source_to_csv/1'
source_file_path_linux = '/share/home/4120107034/MMOSurv/Datasets_and_Preprocessing/data/DNA_methylation_source_to_csv/Linux/1'
target_file_path = '/share/home/4120107034/MMOSurv/Datasets_and_Preprocessing/data/DNA_methylation_variable_filter_10000/1'

flag = True
for root, dirs, files in os.walk(source_file_path):
    for file in files:
        cancer_type = file.split('.')[0].split('-')[1]
        if True:
        # if cancer_type not in cancer_list:
            file_path = str(os.path.join(root, file).encode('utf-8'), 'utf-8')
            temp_class = pd.read_csv(file_path, sep=',')
            temp_class['label'] = file.split('.')[0]
            if flag:
                temp = temp_class
                flag = False
            else:
                temp = temp.append(temp_class)
for root, dirs, files in os.walk(source_file_path_linux):
    for file in files:
        cancer_type = file.split('.')[0].split('-')[1]
        if True:
        # if cancer_type not in cancer_list:
            file_path = str(os.path.join(root, file).encode('utf-8'), 'utf-8')
            temp_class = pd.read_csv(file_path, sep=',')
            temp_class['label'] = file.split('.')[0]
            temp = temp.append(temp_class)
print("DNA_Methylation_feature的数目：" + str(temp.shape[1]))
# temp = temp.dropna(axis=0, how='any')
temp = temp.dropna(axis=1, thresh=temp.shape[0] * 0.8)
print("删除None后" + "DNA_Methylation_feature的数目：" + str(temp.shape[1]))


columns_name = list(temp_class.columns)
var_value_list = []
colname_list = []
for col in columns_name:
    if col not in ['submitter_id', 'label']:
        try:
            temp[col] = temp[col].fillna(temp[col].median())
            var = np.var(np.array(temp[[col]]).squeeze())
            var_value_list.append(var)
            colname_list.append(col)
            with open('DNA_col_name_var' + '.log', 'a') as f:
                f.writelines('col_name:' + col + ',var:' + str(var) + '\n')
        except Exception as e:
            print(e)
ddff = pd.DataFrame({'col_name': colname_list, 'var_value': var_value_list})
var_col_name = np.squeeze(
    ddff.sort_values(by='var_value', ascending=False).iloc[:5000][['col_name']].values).tolist()
print(var_col_name)
# temp['submitter_id'] = temp['submitter_id'].str[0:12]
# temp = temp.drop_duplicates('submitter_id')
var_col_name.append('label')
var_col_name.append('submitter_id')
print("筛选后gene的数目：" + str(len(var_col_name) - 2))
df_DNA_methylation = temp[var_col_name]
print(df_DNA_methylation.shape)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-ACC"])].to_csv(target_file_path + "/TCGA-ACC.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-BLCA"])].to_csv(target_file_path + "/TCGA-BLCA.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-BRCA"])].to_csv(target_file_path + "/TCGA-BRCA.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-CESC"])].to_csv(target_file_path + "/TCGA-CESC.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-CHOL"])].to_csv(target_file_path + "/TCGA-CHOL.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-COAD"])].to_csv(target_file_path + "/TCGA-COAD.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-DLBC"])].to_csv(target_file_path + "/TCGA-DLBC.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-ESCA"])].to_csv(target_file_path + "/TCGA-ESCA.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-GBM"])].to_csv(target_file_path + "/TCGA-GBM.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-HNSC"])].to_csv(target_file_path + "/TCGA-HNSC.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-KICH"])].to_csv(target_file_path + "/TCGA-KICH.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-KIRC"])].to_csv(target_file_path + "/TCGA-KIRC.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-KIRP"])].to_csv(target_file_path + "/TCGA-KIRP.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-LAML"])].to_csv(target_file_path + "/TCGA-LAML.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-LGG"])].to_csv(target_file_path + "/TCGA-LGG.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-LIHC"])].to_csv(target_file_path + "/TCGA-LIHC.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-LUAD"])].to_csv(target_file_path + "/TCGA-LUAD.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-LUSC"])].to_csv(target_file_path + "/TCGA-LUSC.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-MESO"])].to_csv(target_file_path + "/TCGA-MESO.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-OV"])].to_csv(target_file_path + "/TCGA-OV.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-PAAD"])].to_csv(target_file_path + "/TCGA-PAAD.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-PCPG"])].to_csv(target_file_path + "/TCGA-PCPG.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-PRAD"])].to_csv(target_file_path + "/TCGA-PRAD.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-READ"])].to_csv(target_file_path + "/TCGA-READ.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-SARC"])].to_csv(target_file_path + "/TCGA-SARC.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-SKCM"])].to_csv(target_file_path + "/TCGA-SKCM.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-STAD"])].to_csv(target_file_path + "/TCGA-STAD.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-TGCT"])].to_csv(target_file_path + "/TCGA-TGCT.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-THCA"])].to_csv(target_file_path + "/TCGA-THCA.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-THYM"])].to_csv(target_file_path + "/TCGA-THYM.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-UCEC"])].to_csv(target_file_path + "/TCGA-UCEC.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-UCS"])].to_csv(target_file_path + "/TCGA-UCS.csv", index=False, header=True)
df_DNA_methylation[df_DNA_methylation['label'].isin(["TCGA-UVM"])].to_csv(target_file_path + "/TCGA-UVM.csv", index=False, header=True)

# df_DNA_methylation[~(df_DNA_methylation['label'].isin(["TCGA-BRCA"]))].to_csv(target_file_path + "/BRCA_Del.csv", index=False, header=True)
# df_DNA_methylation[~(df_DNA_methylation['label'].isin(["TCGA-MESO"]))].to_csv(target_file_path + "/MESO_Del.csv", index=False, header=True)