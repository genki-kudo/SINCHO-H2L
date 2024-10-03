from script_sincho.extend_score.constant import *
from script_sincho.extend_score.protein_check import *
from script_sincho.extend_score.pqr_and_involveatom import *
from script_sincho.extend_score.term1_hydro_norm import *
from script_sincho.extend_score.term2_maxdist import *
from script_sincho.extend_score.term3_asa import *
from script_sincho.extend_score.druggability_score import *
from script_sincho.extend_score.distance import *
from script_sincho.extend_score.sascore.get_rearranged_smiles import *
from script_sincho.extend_score.sascore.sascore import *
from subprocess import run
import logging

bash=lambda x:run(x,shell=True)

substituent_list =['OC', 'N(C)C', 'O', 'N', 'C', 'C(C)(C)C', 'N(O)O',
                   'CN', 'C(=O)C', 'S(=O)(=O)C', 'C(=O)N', 'S(=O)(=O)N',
                   'CSC', 'c1ccccc1', 'Cc1ccccc1', 'CCc1ccccc1', 'C1CCC1',
                   'C(C1CC1)', 'CC1CCCC1', 'CC1CCCCC1']

def sa_main(extend_idx, pdb_path):
    args = sys.argv
    mol_dir = os.path.dirname(pdb_path)

    match,smi = read_mol(pdb_path)

    obabel_num = get_obabel_num(match,extend_idx)

    smiles = rearrange_smiles(smi,obabel_num)

    t_file('./subst.smi')
    with open('./subst.smi','a')as out:
        print('smiles', file=out)
    for i in substituent_list:
        ext_smi = smiles+i
        with open('./subst.smi','a')as out:
            print(ext_smi, "avoid-warning-str", file=out)
    
    readFragmentScores("fpscores")

    input_smi = './subst.smi'
    out_path = './out_sa.csv'

    df = pd.read_csv(input_smi)
    suppl = Chem.SmilesMolSupplier(input_smi, sanitize=False)
    sascores = processMols(suppl)
    df['sascore'] = sascores
    #for i in df['sascore']:
        #print(i)
    #print(np.mean(df['sascore']))
    sa_ave = np.mean(df['sascore'])

    #df.to_csv(out_path)

    return sa_ave

def sa_base(pdb_path):
    args = sys.argv
    mol_dir = os.path.dirname(pdb_path)

    match,smi = read_mol(pdb_path)

    t_file('./subst.smi')
    with open('./subst.smi','a')as out:
        print('smiles', file=out)
    with open('./subst.smi','a')as out:
        print(smi, "avoid-warning-str", file=out)
    
    readFragmentScores("fpscores")

    input_smi = './subst.smi'
    out_path = './out_sa.csv'

    df = pd.read_csv(input_smi)
    suppl = Chem.SmilesMolSupplier(input_smi)
    sascores = processMols(suppl)
    df['sascore'] = sascores
    #for i in df['sascore']:
        #print(i)
    #print(np.mean(df['sascore']))
    sa_ave = np.mean(df['sascore'])

    #df.to_csv(out_path)

    return sa_ave
