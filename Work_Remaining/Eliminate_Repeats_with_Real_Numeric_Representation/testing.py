import pandas as pd
import numpy as np
from os import getcwd, listdir


loc_genus = getcwd() + "/Work_Remaining/Eliminate_Repeats_with_Real_Numeric_Representation/real_Genus_100.csv"
df_genus = pd.read_csv(loc_genus)
df_genus

loc_realm = getcwd() + "/Work_Remaining/Eliminate_Repeats_with_Real_Numeric_Representation/real_Realm_100.csv"
df_realm = pd.read_csv(loc_realm)
df_realm

varidnaviria_avgmag = df_realm[df_realm["Sublevel Name"] == "Varidnaviria"]["avg_magnitude"]
duplodnaviria_avgmag = df_realm[df_realm["Sublevel Name"] == "Duplodnaviria"]["avg_magnitude"]
monodnaviria_avgmag = df_realm[df_realm["Sublevel Name"] == "Monodnaviria"]["avg_magnitude"]
riboviria_avgmag = df_realm[df_realm["Sublevel Name"] == "Riboviria"]["avg_magnitude"]


varidnaviria_magtropy = df_realm[df_realm["Sublevel Name"] == "Varidnaviria"]["magtropy"]
duplodnaviria_magtropy = df_realm[df_realm["Sublevel Name"] == "Duplodnaviria"]["magtropy"]
monodnaviria_magtropy = df_realm[df_realm["Sublevel Name"] == "Monodnaviria"]["magtropy"]
riboviria_magtropy = df_realm[df_realm["Sublevel Name"] == "Riboviria"]["magtropy"]

#mag avg
varidnaviria_avgmag.describe()
duplodnaviria_avgmag.describe()
riboviria_avgmag.describe()
monodnaviria_avgmag.describe()


#magtropy
varidnaviria_magtropy.describe()
duplodnaviria_magtropy.describe()
monodnaviria_magtropy.describe()
riboviria_magtropy.describe()


df_genus_1 = df_genus.drop(["Sublevel Name"], axis = 1)
df_realm_1 = df_realm.drop(["Sublevel Name"], axis = 1)


df = pd.concat([df_genus_1, df_realm_1])

df = df.reset_index(drop=True)

df_gpby = df.groupby(list(df.columns))

idx = [x[0] for x in df_gpby.groups.values() if len(x) != 1]
df.reindex(idx)



for row in range(len(df_genus_1)):
    print((df_genus_1.iloc[row, ].values == df_realm_1.values).all())
