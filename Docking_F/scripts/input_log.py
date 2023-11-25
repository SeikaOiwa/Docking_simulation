import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='ログファイルの生成（log.csv）',epilog='Module: pandas,argparse')
parser.add_argument('log_save_Fpath', type=str, help='logファイルの保管フォルダパス')
parser.add_argument('col_n', type=str, help='カラム名（"col1 col2 col3"）いくつでもOK')
parser.add_argument('input_data', type=str, help='入力情報（"hoge hoho ddd"）データ点数＝カラム数とすること')

args = parser.parse_args()
log_save_Fpath = args.log_save_Fpath
col_n = [str(i) for i in args.col_n.split(" ")]
input_data = [str(k) for k in args.input_data.split(" ")]

# log生成
df = pd.DataFrame(data=[input_data],columns=[col_n])

# log保存
df.to_csv(f'{log_save_Fpath}/log.csv',index=False)


