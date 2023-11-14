import pandas as pd
import os
import sys

log_file_path = sys.argv[1]
col1_n = sys.argv[2]
col2_n = sys.argv[3]
col3_n = sys.argv[4]
col1_data = sys.argv[5]
col2_data = sys.argv[6]
input_data = sys.argv[7]

log = pd.read_csv(f'{log_file_path}')

for i in range(len(log)):
   if (log.loc[i,col1_n]==col1_data) and (log.loc[i,col2_n]==col2_data):
      log.loc[i,col3_n]= input_data

   log.to_csv(f'{log_file_path}',index=False)

