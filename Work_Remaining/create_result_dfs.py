
import numpy as np
import pandas as pd
from os import getcwd, listdir, system


taxonomic_level = input("taxonomic level: ")
path = getcwd() + f"/Work_Remaining/covid_EIIP_results/{taxonomic_level}"
csv_files = sorted(listdir(path))
csv_files
entropy_df = pd.read_csv(f"{path}/{csv_files[0]}")
entropy_df.insert(loc=0, column='feature', value='entropy')

magnitude_df = pd.read_csv(f"{path}/{csv_files[3]}")
magnitude_df.insert(loc=0, column='feature', value='magnitude')

magtropy_df = pd.read_csv(f"{path}/{csv_files[4]}")
magtropy_df.insert(loc=0, column='feature', value='magtropy')

fd_df = pd.read_csv(f"{path}/{csv_files[2]}")
fd_df.insert(loc=0, column='feature', value='fd')

fd_magtropy_df = pd.read_csv(f"{path}/{csv_files[1]}")
fd_magtropy_df.insert(loc=0, column='feature', value='fd+magtropy')


result_lst = [entropy_df, magnitude_df, magtropy_df, fd_df, fd_magtropy_df]

result_df = pd.concat(result_lst)
result_df = result_df[result_df['tags.Source'] == 'create_model']
result_df = result_df.loc[0]
result_df

saved_path = getcwd() + f"/Work_Remaining/covid_EIIP_results"

result_df.to_csv(saved_path + f"/{taxonomic_level}_feature_metrics.csv")

# pd.read_csv(f"{saved_path}/{taxonomic_level}_feature_metrics.csv")
