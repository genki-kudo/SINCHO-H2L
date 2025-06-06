from script_sincho.basic.basic_func import *
from Bio.PDB import PDBParser, PDBIO, Select, Structure, Model
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP


def shrakerupley_define():
    sr = ShrakeRupley()
    custom_radii_dict = sr.radii_dict.copy()
    custom_radii_dict.update({"LI":1.82, "NE":1.54, "MN":1.61, "CO":1.52, "FE":1.56, 
                              "AR":1.88, "GA":1.87, "AS":1.85, "KR":2.02, "PD":1.63,
                              "AG":1.72, "IN":1.93, "SN":2.17, "TE":2.06, "XE":2.16,
                              "PT":1.75, "AU":1.66, "HG":1.55, "TL":1.96, "PB":2.02,
                              })
    sr = ShrakeRupley(radii_dict=custom_radii_dict)
    #print(custom_radii_dict)
    return sr, custom_radii_dict


def dist_calc(ligand, poc_and_atom, clusterdir):
    poc = clusterdir+poc_and_atom[0]
    atom = poc_and_atom[1].split('_')[-1]

    for i in open(ligand,'r').readlines():
        if i[12:16].replace(' ','')==atom:
            vec_atom = vec_xyz(i)
    pqr_inf = []
    for i in open(poc,'r').readlines():
        pqr_inf.append(vec_xyz(i))
    
    distance_min = 50000000
    for i in pqr_inf:
        distance_min = dist_cf(vec_atom, i, distance_min)
    
    return distance_min


def fix_elements(structure, custom_radii_dict):
    """
    原子の element 属性を get_name から推定して修正する。
    """
    for atom in structure.get_atoms():
        name = atom.get_name().strip()
        if atom.element is None or atom.element.upper() not in custom_radii_dict:
            # 原子名から元素記号を推定（例：HOHの" O " → O）
            if name.startswith(('CL', 'BR', 'NA', 'ZN', 'FE', 'MG', 'CA')):
                atom.element = name[:2].upper()
            else:
                atom.element = name[0].upper()
    return structure

def around_atoms_search(pqrfile, protein, dist, outputdir, ligand):

    #PDBmoduleで周辺原子をサーチ＋記述子計算
    class AtomSelect(Select):
        def __init__(self, atoms):
            self.atoms = atoms
        def accept_atom(self, atom):
            return atom in self.atoms

    output_pdb = outputdir+"/pocket_environment/pocenv"+pqrfile.split("/")[-1].split(".")[0][7:]+"_"+str(dist)+".pdb"

    # 残基の絞り込み
    accept_residues = ["ALA", "GLY", "VAL", "LEU", "ILE", "MET", "PHE", "TYR", "TRP", "PRO",
                       "SER", "THR", "CYS", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS"]
    convert_residues = {"HIE": "HIS", "HID": "HIS", "HIP": "HIS", "CYX": "CYS", "CYM": "CYS",
                        "ASH": "ASP", "GLH": "GLU"}

    parser = PDBParser(QUIET=True)

    #タンパク質の原子取得
    protein_structure = parser.get_structure('protein', protein)
    protein_atoms = list(protein_structure.get_atoms())


    #ポケットの原子取得
    pocket_pdbfile = "sincho-output/"+pqrfile.split("/")[-1]+".pdb"
    pocket = parser.get_structure('pocket', pocket_pdbfile)

    # ポケットくらすたー及びリガンド原子のすべての座標を取得
    pocket_atoms = []
    for i in open(pocket_pdbfile):
        if i[0:6]=="ATOM  " or i[0:6]=="HETATM":
            pocket_atoms.append((ext_xyz(i)))
    #print(pocket_atoms)
    for i in open(ligand):
        if i[0:6]=="ATOM  " or i[0:6]=="HETATM":
            pocket_atoms.append((ext_xyz(i)))


    #化合物の${dist}Å以内の原子を取得ー＞close_atoms_refine
    neighbors_search = NeighborSearch(protein_atoms)
    close_atoms = set()
    for atom in pocket_atoms:
        close_atoms.update(neighbors_search.search(atom, dist))
    #print(close_atoms)
    close_atoms_refine = set()
    close_atoms_refine_ids = []
    for atom in close_atoms:
        if atom.get_parent().get_resname() in convert_residues:
            atom.get_parent().resname = convert_residues[atom.get_parent().get_resname()]
        if atom.get_parent().get_resname() in accept_residues:
            close_atoms_refine.add(atom)
            close_atoms_refine_ids.append(atom.get_serial_number())
    #ポケット原子をPDBファイルとして書き出し
    io = PDBIO()
    io.set_structure(protein_structure)
    io.save(output_pdb, AtomSelect(close_atoms_refine))
    #ポケットのMolオブジェクトとして読み込み
    close_mol = Chem.MolFromPDBFile(output_pdb, removeHs=False)

    sr, custom_radii_dict = shrakerupley_define()
    rds = [a.element for a in list(protein_structure.get_atoms())]
    #print(set(rds))
    restructure = fix_elements(protein_structure, custom_radii_dict )
    sr.compute(restructure, level="A")

    #logP計算
    poc_logP = MolLogP(close_mol)
    pocket_sasa = 0
    pocket_psa = 0
    for atom in protein_structure.get_atoms():
        if atom.get_serial_number() in close_atoms_refine_ids:
            pocket_sasa += atom.sasa
            if atom.element not in ["C","H"]:
                pocket_psa += atom.sasa

    return poc_logP, pocket_sasa, pocket_psa
    


    


