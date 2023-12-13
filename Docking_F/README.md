# Docking_Simulation
autodock_vinaによるドッキングシミュレーション

## 1. 背景
・近年、医薬品開発において、構造ベース創薬（SBDD：Structure-Based Drug Design）が盛んに検討されている
・農薬開発においても、SBDDを用いた活性本体の探索が検討されている
・Docking_simulationはオープンソースである[autodock vina](https://github.com/ccsb-scripps/AutoDock-Vina)をベースに構成している。

## 2. 仕組み
・インプットデータは酵素3Dデータ（.pdb）とLigandデータ(smiles)
・Docking_simulation内で酵素、Ligandデータに水素、電荷を付加
・酵素の位置情報は固定した上でLigandのポージングを様々に変化させた場合の自由Eを計算し、自由Eが最小値を探索
![イメージ](./file/fig.png)

## 3. 環境構築
(1) Conda環境設定

`conda create -n docking python`

`conda activate -n docking`

(2) 基本環境

`conda install -c conda-forge numpy swig boost-cpp sphinx sphinx_rtd_theme xlsxwriter`

`conda install -c conda-forge xlsxwriter`

`pip install pandas`

`pip install streamlit`

`pip install streamlit-aggrid`

`pip install rdkit`

(3) Docking_simulation tool

`pip install vina`

`conda install -c conda-forge openbabel`

`pip install oddt`

`pip install pubchempy`

## 4.フォルダ構造

・Docking_simulationはフォルダ構造に依存するため下図の通りにフォルダ作成する必要あり

・gitの「Docking_S」フォルダをダウンロードして使用を推奨

![folder](./file/Folder_structure.png)

## 5.使用方法

(1) 作業環境構築

gitから「BGC_pred」フォルダをダウンロード

(2) streamlit起動

BGC_Predフォルダ下で、下記コマンドを実行（local_hostが立ち上がる）

`streamlit run st_pred_bgc.py`
