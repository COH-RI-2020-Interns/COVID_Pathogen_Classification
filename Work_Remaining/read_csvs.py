
import numpy as np
import pandas as pd
from os import getcwd, listdir, system

path = getcwd() + f"/fd_magtropy_csvs"

csv_files = listdir(path)
csv_files = sorted(csv_files[:-1])

df = pd.read_csv(f"{path}/{csv_files[0]}")
df
df['Sublevel Name'].value_counts()

taxonomic_level = input("taxonomic level: ")
path2 = getcwd() + f"/Work_Remaining/betacoronavirus_real_results/{taxonomic_level}"
csv_files2 = sorted(listdir(path2))
csv_files2
entropy_df = pd.read_csv(f"{path2}/{csv_files2[0]}")
entropy_df.insert(loc=0, column='feature', value='entropy')
# entropy_result = entropy_df.iloc[0]

magnitude_df = pd.read_csv(f"{path2}/{csv_files2[3]}")
magnitude_df.insert(loc=0, column='feature', value='magnitude')

magtropy_df = pd.read_csv(f"{path2}/{csv_files2[4]}")
magtropy_df.insert(loc=0, column='feature', value='magtropy')

fd_df = pd.read_csv(f"{path2}/{csv_files2[2]}")
fd_df.insert(loc=0, column='feature', value='fd')

fd_magtropy_df = pd.read_csv(f"{path2}/{csv_files2[1]}")
fd_magtropy_df.insert(loc=0, column='feature', value='fd+magtropy')


result_lst = [entropy_df, magnitude_df, magtropy_df, fd_df, fd_magtropy_df]

result_df = pd.concat(result_lst)
result_df = result_df[result_df['tags.Source'] == 'create_model']
result_df

saved_path = getcwd() + f"/Work_Remaining/betacoronavirus_real_results"

result_df.to_csv(saved_path + f"/{taxonomic_level[2:]}_feature_metrics.csv")

pd.read_csv(f"{saved_path}/{taxonomic_level[2:]}_feature_metrics.csv")
