import streamlit as st
from st_aggrid import GridOptionsBuilder, AgGrid, JsCode
from st_aggrid.shared import ColumnsAutoSizeMode
import pandas as pd
import glob
import os
import shutil
import glob
import os
from io import BytesIO
import oddt
from oddt import toolkit
from oddt import docking
from rdkit import Chem
from IPython.display import SVG
from rdkit.Chem import AllChem, Draw, Descriptors, PandasTools
from vina import Vina
import subprocess
import time
import numpy as np
import math

base_path = os.getcwd()
log_path = f'{base_path}//log'
recept_path =f'{base_path}/enzyme'
ligand_path = f'{base_path}/ligand'
file_path = f'{base_path}/scripts'
result_path = f'{base_path}/Result'
fig_path = f'{base_path}/static'
download_path = base_path
figure = f'{base_path}/file'

def ready_enzyme(type):
   st.markdown(f"""
   ## STEP1:酵素データの準備 
   ### 1) 酵素の3D structureファイル({type})を取得
   - [Protein data bank](https://www.rcsb.org/): 実験に基づく3D構造データ
   - [Alpha fold2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb): 機械学習による立体構造予測サイト
   ### 2) 3D strucureファイルを格納（複数可）
   - notebooks/Mydata_analysis/DockingSimulation/receptor
   ### 3) docking_simulation用データの作成
   - :blue[水素]、:orange[電荷]を付加したpdbqtファイルに変換
   """)

   go = st.button('変換開始')
   if go:
      for enzyme in glob.glob(f'{recept_path}/*.{type}'):
         fname = os.path.splitext(os.path.basename(enzyme))[0]
         os.makedirs(f'{recept_path}/{fname}',exist_ok = True)         
         protein = next(oddt.toolkit.readfile(type,enzyme))
         oddt.docking.AutodockVina.write_vina_pdbqt(protein, recept_path, flexible = False)
         os.remove(enzyme)
         
         for rfpath in glob.glob(f'{recept_path}/*.pdbqt'):
            os.rename(rfpath, f'{recept_path}/{fname}/{fname}.pdbqt')
   
   st.markdown("""### 4) データ削除""")
   flist = os.listdir(f'{recept_path}')
   add = ['-選択-']
   fflist = add+flist
   folder = st.selectbox('データ選択',fflist)
   gogo = st.button('削除')
   if gogo:
      shutil.rmtree(f'{recept_path}/{folder}')
                  
def convert_smiles(smile,name):
   mh = Chem.AddHs(Chem.MolFromSmiles(str(smile)))
   AllChem.EmbedMolecule(mh, AllChem.ETKDGv2()) 
   writer = Chem.SDWriter(f'{ligand_path}/{str(name)}.sdf')
   writer.write(mh)   
   
   if os.path.isfile(f'{ligand_path}/{str(name)}.sdf'): 
      ligands = next(oddt.toolkit.readfile('sdf',f'{ligand_path}/{str(name)}.sdf')) 
      oddt.docking.AutodockVina.write_vina_pdbqt(ligands,ligand_path2,name_id=f'{str(name)}')
      os.remove(f'{ligand_path}/{str(name)}.sdf') 
      
def convert_sdf(name,ligand_path2):
   ligands = next(oddt.toolkit.readfile('sdf',f'{ligand_path}/{str(name)}.sdf'))
   oddt.docking.AutodockVina.write_vina_pdbqt(ligands,ligand_path2,name_id=f'{str(name)}')
   os.remove(f'{ligand_path}/{str(name)}.sdf') 
                              
def ready_ligand_smiles():
   st.markdown("""
   ## STEP2: Ligandデータの準備 (:red[Smiles])
   #### 1) Smilesデータの入手
   - [PubChem](https://pubchem.ncbi.nlm.nih.gov/) 
   #### 2) データ入力
   - :green[直接入力]
   """)
                                                           
   col1,col2 = st.columns((5,10))                           
   with col1:
        lname_1 = st.text_input('化合物名','化合物1',key = 'a')
   with col2:          
        smile_1 = st.text_input('データ１',value = '入力',key = 'a2')
        t_data1 = [lname_1,smile_1]
                        
   col3,col4 = st.columns((5,10))                           
   with col3:
        lname_2 = st.text_input('化合物名','化合物2',key = 'b')
   with col4:          
        smile_2 = st.text_input('データ2',value = '入力',key = 'b2') 
        t_data2 = [lname_2,smile_2]
                                     
   col5,col6 = st.columns((5,10))                           
   with col5:
        lname_3 = st.text_input('化合物名','化合物3',key = 'c')
   with col6:          
        smile_3 = st.text_input('データ3',value = '入力',key = 'c2') 
        t_data3 = [lname_3,smile_3]
                       
   col7,col8 = st.columns((5,10))                           
   with col7:
        lname_4 = st.text_input('化合物名','化合物4',key = 'd')
   with col8:          
        smile_4 = st.text_input('データ4',value = '入力',key = 'd2')
        t_data4 = [lname_4,smile_4]
            
   col9,col10 = st.columns((5,10))                           
   with col9:
        lname_5 = st.text_input('化合物名','化合物5',key = 'e')
   with col10:          
        smile_5 = st.text_input('データ5',value = '入力',key = 'e2')
        t_data5 = [lname_5,smile_5]
     
     
   st.markdown(""" - :green[アップロード] """)
      
   df = pd.read_csv(f'{figure}/smiles.csv')
   df.to_csv(buf := BytesIO(), index=True)
   st.download_button("雛型ダウンロード",buf.getvalue(),'smiles.csv',"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
   smile_list = st.file_uploader('smileデータ',type='csv',key='f')
   if smile_list:
      df2 = pd.read_csv(smile_list)
      st.dataframe(df2)
   
   st.markdown(""" #### 3) ligandデータ保存名の設定""")
   exit_list = os.listdir(f'{ligand_path}')
   lf_name = st.text_input('')
   if lf_name in exit_list:
      st.markdown("""#### :red[既に存在します。名称変更してください]""")
   
   st.markdown("""
   #### 4)docking_simulation用データの作成
   - :blue[水素]、:orange[電荷]を付加したpdbqtファイルに変換 
   """)
   
   go2 = st.button('変換開始')
   if go2:
 
      ligand_path2 = f'{ligand_path}/{lf_name}'
      os.makedirs(ligand_path2,exist_ok = True)
      smiles_data_list = []
    
      if smile_1 != '入力':
        convert_smiles(smile_1,lname_1,ligand_path2)
        smiles_data_list.append(t_data1)
      if smile_2 != '入力':
        convert_smiles(smile_2,lname_2,ligand_path2)
        smiles_data_list.append(t_data2)
      if smile_3 != '入力':
        convert_smiles(smile_3,lname_3,ligand_path2)
        smiles_data_list.append(t_data3)
      if smile_4 != '入力':
        convert_smiles(smile_4,lname_4,ligand_path2)
        smiles_data_list.append(t_data4)
      if smile_5 != '入力':
        convert_smiles(smile_5,lname_5,ligand_path2)
        smiles_data_list.append(t_data5)
      
      if not smiles_data_list:
         print('no_data')
      else:
         smiles_data = pd.DataFrame(smiles_data_list,columns=[['name','smiles']])
         smiles_data.to_csv(f'{ligand_path2}/smile_list.csv')
       
      if smile_list:
         df2.to_csv(f'{ligand_path2}/smile_list.csv')         
         for i in range(0,len(df2)):
            s_name = df2.loc[i,'name']
            smiles_data = df2.loc[i,'smiles']
            subprocess.run(["python",f"{file_path}/convert_sdf.py",s_name,smiles_data,ligand_path2])            
            subprocess.run(["python",f"{file_path}/convert_pdbqt.py",ligand_path2])
            #convert_smiles(smiles_data,s_name)
            
   st.markdown(""" #### 5)データ削除 """)
   l_list = os.listdir(f'{ligand_path}')
   select_l = st.selectbox('削除対象を選択',l_list)
   elase = st.button('削除')
   if elase:
      ligand_path3 = f'{ligand_path}/{select_l}'
      shutil.rmtree(ligand_path3)           

def ready_ligand_sdf():
   st.markdown("""
   ## STEP2: Ligandデータの準備 (:red[sdfファイル])
   #### 1) sdfデータの入手
   - [PubChem](https://pubchem.ncbi.nlm.nih.gov/) 
   #### 2) データ保存(複数可)
   - notebooks/DockingSimulation/ligand
   #### 3) ligandデータ保存名の設定""")
  
   exit_list = os.listdir(f'{ligand_path}')
   lf_name = st.text_input('')
   if lf_name in exit_list:
      st.markdown("""#### :red[既に存在します。名称変更してください]""")
     
   st.markdown("""
   #### 4)docking_simulation用データの作成
   - :blue[水素]、:orange[電荷]を付加したpdbqtファイルに変換  
   """)
  
   go3 = st.button('変換開始')
   if go3:
      for i in glob.glob(f'{ligand_path}/*.sdf'):
         name = os.path.splitext(os.path.basename(i))[0]
         convert_sdf(name,lf_name)  
   
   st.markdown(""" #### 5)データ削除 """)
   l_list = os.listdir(f'{ligand_path}')
   select_l = st.selectbox('削除対象を選択',l_list)
   elase = st.button('削除')
   if elase:
      ligand_path3 = f'{ligand_path}/{select_l}'
      shutil.rmtree(ligand_path3)

def docking_condition():
   st.markdown("""
   ## STEP3: Docking_simulation （:red[条件設定]） 
   #### 酵素毎に条件設定(.CSV)を作成
   - 1) 基質結合部位の中心座標（予測値）[ProteinsPlus](https://proteins.plus/)
   - 2) 探索範囲（x,y,z）
   - 3) 探索回数(exhaustiveness) ※低いと予測精度が落ちる（50>に設定)
   - 4) 取得データ数 10コ程度
   """)
   
   df = pd.read_csv(f'{figure}/docking_condition.csv')
   df.to_csv(buf := BytesIO(), index=True)
   st.download_button("雛型ダウンロード",buf.getvalue(),'docking_condition.csv',"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
   st.image(f'{figure}/fig2.png')
   st.markdown("""#### アップロード""")
   condition_data = st.file_uploader('アップローダー',type='csv',key='f')
   if condition_data:
      df = pd.read_csv(condition_data)
      st.dataframe(df)
      st.markdown("""#### 酵素の選択 """)
      flist = os.listdir(f'{recept_path}')
      add = ['-選択-']
      fflist = add+flist
      folder = st.selectbox('データ選択',fflist)
      gogo = st.button('保存')
      if gogo:
         df.to_csv(f'{recept_path}/{folder}/docking_condition.csv')

    
def docking_simulation():
   st.markdown("""
   ## STEP3: Docking_simulation （:red[解析]）
   - 出力ファイル名の入力
   """)
   
   fname = st.text_input('(任意)','docking_result')
   
   process_num = st.slider('プロセス数',1,20,1,1)
   max_search_time_ = st.slider('解析時間の上限(分)',1,60,20,1)
   max_search_time = max_search_time_ * 60

   st.markdown(""" - 酵素の選択 """)
   flist = os.listdir(f'{recept_path}')
   enzyme_name = st.multiselect('選択',flist)
   
   st.markdown(""" - ligandリストの選択 """)
   add = ['-選択-']
   l_list_p = os.listdir(f'{ligand_path}')
   l_list = add + l_list_p
   select_ligand = st.selectbox('選択',l_list)
   select_ligand_path = f'{ligand_path}/{select_ligand}'
   if select_ligand != '-選択-':
      llist = os.listdir(select_ligand_path)
      st.dataframe(llist)
      
   st.markdown(""" - シミュレーション開始 """)
   go_simulation = st.button('開始')

   # docking
   if go_simulation:
      # 進捗確認用
      os.makedirs(f'{result_path}/{fname}',exist_ok=True)      
      en_num = len(enzyme_name) 
      lig_num = len([i for i in glob.glob(f'{select_ligand_path}/*.pdbqt',recursive=True)])
      total_analysis_num = en_num * lig_num
      with open(f'{result_path}/{fname}/total_analysis_num.txt','a') as f:
         f.write(f'total_analysis_num {total_analysis_num}')

      for en in enzyme_name:
         os.makedirs(f'{result_path}/{fname}/{en}',exist_ok=True)
         shutil.copy2(f'{select_ligand_path}/smile_list.csv',f'{result_path}/{fname}/smile_list.csv')       
         for lig_path in glob.glob(f'{select_ligand_path}/*.pdbqt'): 
            lig_name = os.path.splitext(os.path.basename(lig_path))[0]  

            # docking_引数
            denzyme_path = f'{recept_path}/{en}/{en}.pdbqt'
            doc_path =  f'{recept_path}/{en}/docking_condition.csv'
            dligand_path = lig_path
            save_path = f'{result_path}/{fname}'

            # docking_simulation(multi process)
            procs = {}
            while sum(proc.poll() is None for n,proc in procs.items())>=int(process_num):
               time.sleep(10)
            proc = subprocess.Popen(["python",f"{file_path}/pre_docking.py",
            en,
            lig_name,
            denzyme_path,
            doc_path,
            dligand_path,
            save_path,
            file_path,
            "--max_time",str(max_search_time)
            ],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

            procs[f'{en}_{lig_name}']= proc

def confirm_progress():
   st.markdown("""
   ## STEP3: 解析状況の確認
   - Docking simulationデータの選択
   """)
   result_list = os.listdir(f'{result_path}')
   add = ['-選択-']
   result_list2 = add + result_list
   select_result = st.selectbox('解析データ選択',result_list2)

   go = st.button('解析状況の確認')
   if go:
      # 解析総数の確認
      with open(f'{result_path}/{select_result}/total_analysis_num.txt','r') as f:
         num = f.read()
         total_num = int(num.split(" ")[1])

      # 進捗状況まとめ
      # 完了
      success_num = len([i for i in glob.glob(f'{result_path}/{select_result}/**/*/*.pdbqt',recursive=True)])
      # time_out or error
      logs = [pd.read_csv(log) for log in glob.glob(f'{result_path}/**/*/*/log.csv',recursive=True)]
      log_df_ = pd.concat(logs,axis=0)
      log_df = log_df_.reset_index()
      timeout_num = len(log_df[log_df['Progress'].str.contains('Time_out')])
      error_num = len(log_df[log_df['Progress'].str.contains('Error')])

      finish_ratio = math.floor((success_num + timeout_num +error_num)/total_num *100)

      # 進捗状況表示
      st.write(f'進捗率：{finish_ratio} %')
      st.progress(finish_ratio,text='')
      st.markdown("""#### (内訳)""")
      st.write(f'総解析：{total_num} 個')
      st.write(f'データ生成：{success_num} 個 ({round(success_num/total_num*100,1)}%)')
      st.write(f'エラー：{error_num} 個 ({round(error_num/total_num*100,1)}%)')
      st.write(f'解析中断：{timeout_num} 個 ({round(timeout_num/total_num*100,1)}%)')
         
def view_structure(dataframe):
 
   #画像保存フォルダの更新
   shutil.rmtree(f'{fig_path}')
   os.makedirs(f'{fig_path}')
 
   describe_fig_path =[] 
   
   for i in range(len(dataframe)):
      #要素抜出し
      smile = dataframe.loc[i,'smiles']
      lig_name = dataframe.loc[i,'ligand']
      mol = Chem.MolFromSmiles(smile)
      
      #保存先のパス
      fig_path2 = f'{fig_path}/{lig_name}.png'
      
      #画像表示用のパス
      path = f'./app/static/{lig_name}.png'
      describe_fig_path.append(path)
      
      #smile⇒二次構造
      Draw.MolToFile(mol,fig_path2,size=(100, 100))
      
   #データフレームにパス追加
   dataframe2 = dataframe[['Enzyme','ligand','dG']]
   dataframe2['apps'] = describe_fig_path   
   #st.dataframe(dataframe2)
   
   #画像入りデータフレームの表示
   st.data_editor(
    dataframe2,
    column_config={"apps": st.column_config.ImageColumn("Preview Image", help="Streamlit app preview screenshots")},
    hide_index=True,)
      
def calculate_data(select_result,select_enzyme):
   
   all_docking = []
     
   for i in glob.glob(f'{result_path}/{select_result}/{select_enzyme}/**/*.pdbqt',recursive=True):         
      df = pd.read_table(i)
      df2 = df[df['MODEL 1'].str.contains('RESULT')]
      df3 = df2.reset_index(drop=True)
      df4 = df3['MODEL 1'].str.split(':',expand=True)
      df5 = df4[1].str.split('      ',expand=True)

      ligand_name = os.path.splitext(os.path.basename(i))[0]    
      min_delta = df5.loc[0,0]
      add_data = [select_enzyme,ligand_name,min_delta]
      all_docking.append(add_data)

   result = pd.DataFrame(all_docking)
   result = result.set_axis(['Enzyme','ligand','dG'],axis='columns',copy=False).sort_values('dG',ascending=False)
     
   add_smile = pd.read_csv(f'{result_path}/{select_result}/smile_list.csv')
   
   merge_data = pd.merge(result,add_smile, left_on='ligand',right_on='name', how='left')
   merge_data2 = merge_data[['Enzyme','ligand','dG','smiles']]
   merge_data3 = merge_data[['Enzyme','ligand','dG','smiles']]
     
   #data表示
   st.dataframe(merge_data2)
   
   #ダウンロードデータの生成
   go = st.button('ダウンロード準備')
   if go:
   
      #smilesから構造データ構築
      PandasTools.AddMoleculeColumnToFrame(merge_data2, molCol='Structure', smilesCol='smiles')
   
      #excelファイル保存   
      PandasTools.SaveXlsxFromFrame(merge_data2, f"{download_path}/{select_enzyme}_ligand_Docking.xlsx", 
                                 molCol='Structure', 
                                 size=(150,150)
                                 )
      #streamlit上で結果概要を表示
      view_structure(merge_data3)
     
      st.markdown(""" #### 下記リンクをコピーしブラウザに貼付けてください """)
      st.code(f"file:///D:/Box/[Azure_box]_MI-ATRL-3/bio_lab/Docking_simulation/{select_enzyme}_ligand_Docking.xlsx")
      
      merge_data3.to_excel(buf := BytesIO(), index=True)
      st.download_button("Download (画像なし)",buf.getvalue(),f"{select_enzyme}_ligand_Docking.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
         
def confirm_result():
   st.markdown("""
   ## STEP4: 解析結果の確認
   - Docking simulationデータの選択
   """)
   result_list = os.listdir(f'{result_path}')
   add = ['-選択-']
   result_list2 = add + result_list
   select_result = st.selectbox('解析データ選択',result_list2)
   
   if select_result != '-選択-':
      result_enzyme_list_p = os.listdir(f'{result_path}/{select_result}')
      result_enzyme_list_p.remove('smile_list.csv')
      result_enzyme_list = add + result_enzyme_list_p
      select_enzyme = st.selectbox('酵素選択',result_enzyme_list)
      if select_enzyme != '-選択-': 
         calculate_data(select_result,select_enzyme)
        
      elase = st.button('解析データの削除',key='elase')
      if elase:
         shutil.rmtree(f'{result_path}/{select_result}')
         shutil.rmtree(f"{download_path}/{select_enzyme}_ligand_Docking.xlsx")
        

st.sidebar.subheader('STEP1:酵素データの準備')
select = st.sidebar.selectbox('入力データの選択',['-select-','.mol2','.pdb'])
st.sidebar.subheader('STEP2: Ligandデータの準備')
select2 = st.sidebar.selectbox('入力データの選択',['-select-','smiles','.sdf'])
st.sidebar.subheader('STEP3: Docking_simulation')
select3 = st.sidebar.selectbox('選択',['-select-','条件設定','解析','進捗確認'])
st.sidebar.subheader('STEP4: 解析結果の確認')
select4 = st.sidebar.selectbox('選択',['-select-','解析結果の確認'])

st.markdown(""" # :green[Docking Simulation Tool] """)
st.image(f'{figure}/fig.png')

if select == '.mol2':
   ready_enzyme('mol2')
  
if select == '.pdb': 
   ready_enzyme('pdb')
  
if select2 == 'smiles':
   ready_ligand_smiles()
  
if select2 == '.sdf':
   ready_ligand_sdf()
  
if select3 == '条件設定':
   docking_condition()

if select3 == '解析':
   docking_simulation()
   
if select3 == '進捗確認':
   confirm_progress()

if select4 == '解析結果の確認':
   confirm_result()





