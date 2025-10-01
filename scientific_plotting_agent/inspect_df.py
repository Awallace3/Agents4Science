import pandas as pd
import sys
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)
sys.stdout.reconfigure(encoding='utf-8')
df = pd.read_pickle('combined_df_4569.pkl')
print(df.columns.tolist())
