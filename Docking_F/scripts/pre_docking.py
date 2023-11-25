import os
import shutil
import pandas as pd
import argparse
import subprocess
import time

def make_result_f(r_path,en_name,lig_name):
    """docking結果を格納するフォルダ生成、フォルダパス情報を返す
       フォルダ構造：親フォルダ ー 酵素名フォルダ ー リガンド名フォルダ
       酵素×複数リガンドのdockingの場合、同一酵素名の下に複数のリガンド名フォルダが生成 
    Parameters:
    -----------
    r_path:str
        親フォルダのパス
    en_name:str
        酵素フォルダ名
    lig_name:str
        リガンドフォルダ名
    
    Returns:
    --------
    folder_path: str
    """
    folder_path = f'{r_path}/{en_name}/{lig_name}'
    os.makedirs(folder_path,exist_ok=True)

    return folder_path

def trans_docking_condition(dc_path,save_f_path):
    """docking conditionが記載されたcsvファイルを指定先にコピー

    Parameters:
    ----------
    dc_path: str
        コピー元のファイルパス
    save_f_path: str
        コピー先のフォルダパス
    """
    shutil.copy(dc_path,save_f_path)

def read_docking_condition(dcf_path):
    """docking_condition.csvからドッキング条件を抽出

    Parameters:
    ----------
    dcf_path: str
        docking_condition.csvを格納したフォルダパス
    Returns:
    --------
    binding_cite: str
    search_area: str
    exhaustiveness: str
    """
    df = pd.read_csv(f'{dcf_path}/docking_condition.csv',skiprows=1)
    lx = df.loc[0,'lx']
    ly = df.loc[0,'ly']
    lz = df.loc[0,'lz']
    binding_cite = f"{lx} {ly} {lz}"
    sx = df.loc[0,'sx']
    sy = df.loc[0,'sy']
    sz = df.loc[0,'sz']
    search_area = f"{sx} {sy} {sz}"
    exhaustiveness_ = df.loc[0,'ex_num']
    exhaustiveness = str(exhaustiveness_)

    return binding_cite,search_area,exhaustiveness

def simulation(file_path,denzyme_path,dligand_path,binding_cite,search_area,exhaustiveness,data_num,save_path,max_search_time):
    """autdock_vina を使ったdocking simulationを行い、dGを含むtmp.pdbqtが生成
       1回のドッキングシミュレーションの制限時間はsearch_timeで規定、それ以上かかる場合は結合しないものと判断する

    Parameters:
    ----------
    denzyme_path: str
      enzyme data　path (/**.pdbqt)
    dligand_path: str
      ligand data　path (/**.pdbqt)
    binding_cite: list,float
      基質結合位置情報, [x,y,z]
    search_area: list,float
      リガンドの結合の探索範囲 [x,y,z]
    exhaustiveness: int
      探索回数, 低いと予測精度が落ちる（50>に設定)
    data_num: int
      取得データ数
    save_path: str:
      simulation結果の保存のパス（/usr/*/*/***.pdbpt)
    max_search_time: int/str
      Docking simulation (subprocess.run) の実行時間の上限値

    See Also:
    --------
    add_docking_data : Antismash解析データから構築した既知二次代謝物をDocking simulation、結果を追加したデータフレームを返す
    """

    s_time = time.time()
    try:
      re = subprocess.run(["python",f"{file_path}/docking.py",
      denzyme_path,
      dligand_path,
      save_path,
      binding_cite,
      search_area,
      str(exhaustiveness),
      str(data_num)],
      timeout=int(max_search_time),
      capture_output=True, check=True, encoding='utf-8')
        
      e_time = time.time()
      analysis_time = round((e_time-s_time)/60,1)
        
      return 'finish', analysis_time
       
    except subprocess.CalledProcessError as f:
      return 'Error',f.stderr
    
    except subprocess.TimeoutExpired as e:
      return 'time_out',''


parser = argparse.ArgumentParser(description='Docking simulation (Autodock vina)')
parser.add_argument('enzyme',type=str, help='酵素名')
parser.add_argument('ligand',type=str, help='リガンド名')
parser.add_argument('e_path', type=str, help='酵素ファイルのパス(/**.pdbqt)')
parser.add_argument('dc_path', type=str, help='Docking_condition.csvのパス(/docking_condition.csv)')
parser.add_argument('l_path', type=str, help='Ligandファイルのパス(/**.pdbqt)')
parser.add_argument('r_path',type=str, help='結果の保存先フォルダパス')
parser.add_argument('script_path',type=str, help='docking.pyおよびinput_log.pyを保存したフォルダパス')
parser.add_argument('--max_time',type=int, help='解析時間の上限値(秒),デフォルト20分',default=1200)

args = parser.parse_args()
en_name = args.enzyme
lig_name = args.ligand
denzyme_path = args.e_path
dc_path = args.dc_path
dligand_path = args.l_path
r_path = args.r_path
file_path = args.script_path


# 結果格納フォルダ生成
save_f_path = make_result_f(r_path,en_name,lig_name)

# Docking conditionをコピー
trans_docking_condition(dc_path,save_f_path)

# Docking conditionの抽出
binding_cite,search_area,exhaustiveness = read_docking_condition(save_f_path)

# Docking simulationの実行
data_num = 10
save_path = f'{save_f_path}/{lig_name}.pdbqt'
max_search_time = 20*60
result,analysis_time = simulation(file_path,denzyme_path,dligand_path,binding_cite,search_area,exhaustiveness,data_num,save_path,max_search_time)

# log出力
col_n = "Enzyme Ligand Progress"

if result == 'finish':
    subprocess.run(["python",f"{file_path}/input_log.py",
    save_f_path,
    col_n,
    f"{en_name} {lig_name} {str(analysis_time)}"]) 

if result == 'Error':
    subprocess.run(["python",f"{file_path}/input_log.py",
    save_f_path,
    col_n,
    f"{en_name} {lig_name} {str(analysis_time)}"])           

if result == 'time_out':
    subprocess.run(["python",f"{file_path}/input_log.py",
    save_f_path,
    col_n,
    f"{en_name} {lig_name} Time_out"])    
    