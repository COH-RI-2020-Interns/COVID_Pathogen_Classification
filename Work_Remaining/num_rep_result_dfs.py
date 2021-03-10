
import numpy as np
import pandas as pd
from os import getcwd, listdir, system


taxonomic_level = input("folder name: ")
path = getcwd() + f"/Work_Remaining/num_rep_results/{taxonomic_level}"
csv_files = sorted(listdir(path))
csv_files
int1_df = pd.read_csv(f"{path}/{csv_files[1]}")
int1_df.insert(loc=0, column='feature', value='Int1')
int1_df = int1_df[int1_df['tags.Source'] == 'create_model']

int2_df = pd.read_csv(f"{path}/{csv_files[2]}")
int2_df.insert(loc=0, column='feature', value='Int2')
int2_df = int2_df[int2_df['tags.Source'] == 'create_model']

real_df = pd.read_csv(f"{path}/{csv_files[9]}")
real_df.insert(loc=0, column='feature', value='Real')
real_df = real_df[real_df['tags.Source'] == 'create_model']

EIIP_df = pd.read_csv(f"{path}/{csv_files[0]}")
EIIP_df.insert(loc=0, column='feature', value='EIIP')
EIIP_df = EIIP_df[EIIP_df['tags.Source'] == 'create_model']

pp_df = pd.read_csv(f"{path}/{csv_files[8]}")
pp_df.insert(loc=0, column='feature', value='PP')
pp_df = pp_df[pp_df['tags.Source'] == 'create_model']

pn_df = pd.read_csv(f"{path}/{csv_files[7]}")
pn_df.insert(loc=0, column='feature', value='PN')
pn_df = pn_df[pn_df['tags.Source'] == 'create_model']

justA_df = pd.read_csv(f"{path}/{csv_files[3]}")
justA_df.insert(loc=0, column='feature', value='Just A')
justA_df = justA_df[justA_df['tags.Source'] == 'create_model']

justC_df = pd.read_csv(f"{path}/{csv_files[4]}")
justC_df.insert(loc=0, column='feature', value='Just C')
justC_df = justC_df[justC_df['tags.Source'] == 'create_model']

justG_df = pd.read_csv(f"{path}/{csv_files[5]}")
justG_df.insert(loc=0, column='feature', value='Just G')
justG_df = justG_df[justG_df['tags.Source'] == 'create_model']

justT_df = pd.read_csv(f"{path}/{csv_files[6]}")
justT_df.insert(loc=0, column='feature', value='Just T')
justT_df = justT_df[justT_df['tags.Source'] == 'create_model']

result_lst = [int1_df, int2_df, real_df, EIIP_df, pp_df, pn_df, justA_df, justC_df, justG_df, justT_df]

result_df = pd.concat(result_lst)
result_df = result_df.loc[0]
result_df

saved_path = getcwd() + f"/Work_Remaining/num_rep_results"

result_df.to_csv(saved_path + f"/{taxonomic_level}_metrics.csv")
