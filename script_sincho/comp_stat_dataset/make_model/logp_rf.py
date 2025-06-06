import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors, AllChem
import plotly.graph_objects as go
import numpy as np
from sklearn import datasets
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import seaborn as sns
import pandas as pd
from rdkit.Chem.Crippen import MolLogP
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
from sklearn.linear_model import ElasticNet
from scipy.stats import uniform
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.model_selection import RepeatedKFold, cross_val_score
from sklearn.metrics import mean_squared_error, mean_absolute_error, mean_absolute_percentage_error
from scipy.stats import randint
import pickle
import json

"""
def pass_descriptor(mol):
    conformer = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    dist, atom1, atom2 =0, -1, -1
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            pos1 = np.array(conformer.GetAtomPosition(i))
            pos2 = np.array(conformer.GetAtomPosition(j))
            dist_tmp = np.linalg.norm(pos1 - pos2)
            if dist_tmp>=dist:
                dist = dist_tmp
                atom1, atom2 = i,j
    spath = Chem.GetShortestPath(mol, atom1, atom2)
    atoms = mol.GetAtoms()
    pass_mw = 0
    sbl_list = []
    for a in atoms:
        if a.GetIdx() in spath[1:-1]:
            sbl_list.append(a.GetSymbol())
            pass_mw+=a.GetMass()
    return dist, pass_mw

df = pd.read_csv("descriptors3.csv")

names = df["PDB"]
vol = []
mw = []
dist = []
pass_mw = []
for pdb in names:
    print(pdb)
    sdf = "../20241202_PDBbind_v2020_refined/refined-set/"+pdb+"/"+pdb+"_ligand.sdf"
    mol = Chem.SDMolSupplier(sdf,removeHs=False, sanitize=False)[0]
    vol.append(AllChem.ComputeMolVolume(mol, gridSpacing=1.0))
    mw.append(Chem.Descriptors.MolWt(mol))
    d, pm = pass_descriptor(mol)
    dist.append(d)
    pass_mw.append(pm)

df["volume"] = vol
df["mw"] = mw
df["dist"] = dist
df["pass_mw"] = pass_mw

df.to_csv("descriptors4.csv")
"""



df = pd.read_csv("descriptors4.csv")
print("default_dataset", df.shape)
df = df.dropna()

#filiter ON/OFF

for index, line in enumerate(df.iterrows()):
    id = line[1][4]
    #print(id)
    mol = Chem.SDMolSupplier("../20241202_PDBbind_v2020_refined/refined-set/"+id+"/"+id+"_ligand.sdf",removeHs=False, sanitize=False)[0]

    if Chem.rdMolDescriptors.CalcNumLipinskiHBA(mol)>=10 or Chem.rdMolDescriptors.CalcNumLipinskiHBD(mol)>=5 or Chem.Descriptors.MolWt(mol)>=500:
        #or float(line[1][22])/float(line[1][23])>=0.5:
        #print("del",id, float(line[1][22])/float(line[1][23]))
        df.drop( df [ df["PDB"].str.contains(id) ] .index, inplace=True )

y = df["clogP"]
#dfの4列目以降を重回帰のための説明変数として格納
x = df[ ['plogP_4.5', 'plogP_5.0', 'plogP_5.5', 'plogP_6.0', 'plogP_6.5', 'plogP_7.0', 
         'psasa_4.5', 'psasa_5.0', 'psasa_5.5', 'psasa_6.0', 'psasa_6.5', 'psasa_7.0', 
         'ppsa_4.5', 'ppsa_5.0', 'ppsa_5.5', 'ppsa_6.0', 'ppsa_6.5', 'ppsa_7.0'] ]
print(x.shape)
print(x.columns)


# トレーニング・テストデータに分割（再現性のため random_state=42）
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

# ランダム探索のパラメータ空間（広めにとれる）
param_dist = {
    'n_estimators': randint(50, 1000, 50),
    'max_depth': [None] + list(range(10, 101, 10)),
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'max_features': ['auto', 'sqrt']
}

# RandomizedSearchCV の設定
random_search = RandomizedSearchCV(
    estimator=RandomForestRegressor(random_state=42),
    param_distributions=param_dist,
    n_iter=500,  # 試行回数（多いほど精度↑）
    cv=5,
    scoring='r2',
    n_jobs=-1,
    verbose=1,
    random_state=42
)

random_search.fit(x_train, y_train)

# 最適パラメータ & モデル
best_model = random_search.best_estimator_
print("✅ 最適パラメータ（ランダム探索）:")
print(random_search.best_params_)

cv_scores = cross_val_score(best_model, x, y, cv=5, scoring='r2', n_jobs=-1)
print("✅ クロスバリデーション平均R²:", np.mean(cv_scores))

# --- テストデータでの評価 ---
y_pred = best_model.predict(x_test)
r2 = r2_score(y_test, y_pred)
mae = mean_absolute_error(y_test, y_pred)
rmse = mean_squared_error(y_test, y_pred, squared=False)
mape = mean_absolute_percentage_error(y_test, y_pred)

print(f'\n📊 評価指標:')
print(f'R²: {r2:.3f}')
print(f'MAE: {mae:.3f}')
print(f'RMSE: {rmse:.3f}')
print(f'MAPE: {mape:.3%}')

# --- 特徴量重要度 ---
importances = pd.Series(best_model.feature_importances_, index=x.columns)
importances.sort_values(ascending=True).plot(kind='barh', figsize=(8, 6))
plt.title("Feature Importance")
plt.tight_layout()
plt.savefig("feature_importances.png")
plt.show()

final_model = RandomForestRegressor(**random_search.best_params_, random_state=42)
final_model.fit(x,y)

import joblib
joblib.dump(final_model,"logp_rf.pkl")

