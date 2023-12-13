from vina import Vina
import os
import argparse

parser = argparse.ArgumentParser(description='Docking simulation (Autodock vina)')
parser.add_argument('e_path', type=str, help='酵素ファイルのパス(/**.pdbqt)')
parser.add_argument('l_path', type=str, help='Ligandファイルのパス(/**.pdbqt)')
parser.add_argument('o_path',type=str, help='結果の保存先パス(/**.pdbqt)')
parser.add_argument('b', type=str,help='結合座標：x y z')
parser.add_argument('s', type=str,help='探索座標：x y z')
parser.add_argument('e', type=int, help='探索回数',default=50)
parser.add_argument('d',type=int, help='データサンプリング数',default=10)

args = parser.parse_args()

denzyme_path = args.e_path
dligand_path = args.l_path
binding_cite = [float(i) for i in args.b.split(" ")]
search_area = [float(k) for k in args.s.split(" ")]
exhaustiveness = args.e
data_num = args.d
opath = args.o_path

v = Vina(sf_name='vina')
v.set_receptor(denzyme_path)
v.set_ligand_from_file(dligand_path)
v.compute_vina_maps(center=binding_cite, box_size=search_area)
v.dock(exhaustiveness=exhaustiveness,n_poses=data_num)
v.write_poses(f'{opath}',n_poses=data_num,overwrite=True)