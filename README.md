# SINCHO(Site Idenrification and Next CHOice) protocol for Hit-to-Lead
![図1](https://github.com/genki-kudo/SINCHO-H2L/blob/main/sincho_graphics.png)
SINCHO protocol is a method for **the prediction/suggestion of the desirable anchor atom and growth site pair** for the modification of the hit compound in hit-to-lead process.

## Requirements
* **Pocket to Concavity (P2C)**

  P2C is used for site identification.
  [P2C installation is refferred to the repository.](https://github.com/genki-kudo/Pocket-to-Concavity)
  
  (This guide in the P2C repository includes to install P2C, Fpocket2 and Python3. )

* **Python**

  We tested on the python3.7.10.
  It requires to install several modules. Most modules are already installed when P2C installation.
  * numpy(1.21.5)
  * pandas(1.3.5)
  * scikit-learn(1.0.2)
  * scipy(1.7.0)
  * pymol(open-source)
  * Biopython(1.79)
  You should add Biopython module for the SINCHO protocol.
  ```
  conda install Biopython==1.79
  ```


## Installation
* **SINCHO**  
  Download this source code, and set PATH in this directory.  
  ```
  git clone https://github.com/genki-kudo/SINCHO-H2L.git
  cd SINCHO-H2L/
  echo "export PATH=\$PATH:`pwd`/bin" >> ~/.bashrc
  source ~/.bashrc
  ```

## Preparation of Input Files

In the SINCHO protocol, a **protein 3D structure file** and a **ligand 3D structure file** must be prepared as **PDB format**.

* **Protein 3D structure file**

  Element symbol must be included in column 77-78 of each "ATOM" line.
  If the PDB file contains substrates (such as DNA, RNA, and ligands), it is recommended that they be removed so site identification can function properly.
* **Ligand 3D structure file**
  
  Coordinates in the file must be the binding state in the protein 3D structure file. This file name need to correspond to residue name in the PDB file.
  Hydrogen atoms must be added, and element symbol must be included in column 77-78 of each "ATOM" line. "CONECT" lines are needed to include in the file.
  Atom names in the PDB file need to avoid duplication. 

## Running
Two commands, P2C and SINCHO, are required in the execution.

### 1. P2C 
Option details of P2C is shown in [P2C repository](https://github.com/genki-kudo/Pocket-to-Concavity).
**For the execution SINCHO protocol, P2C is performed by the following command:**
```
$ p2c -m LB -p ${protein}.pdb -l ${ligand}.pdb -d 10
```
```${protein}.pdb``` and ```${ligand}.pdb``` are specified to the input files of your interest.

### 2. SINCHO 
To view options of SINCHO.
```
usage: sincho [-h] [-c CLUSTERDIR] [-o OUTPUTDIR] [-f FPOCKETDIR] -l LIGAND -p
              PROTEIN [-log LOGFILENAME] [-n NUMBER_CAND]

optional arguments:
  -h, --help            show this help message and exit
  -c CLUSTERDIR, --clusterdir CLUSTERDIR
                        specify p2c output directory with pocket clusters if
                        neccesary (default: ./p2c_output/cluster/)
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        specify output directory if neccesary (default:
                        ./sincho-output)
  -f FPOCKETDIR, --fpocketdir FPOCKETDIR
                        specify fpocket output directory if neccesary
                        (default: ./asphere_output)
  -l LIGAND, --ligand LIGAND
                        specify ligand pdb file.
  -p PROTEIN, --protein PROTEIN
                        specify protein pdb file.
  -log LOGFILENAME, --logfilename LOGFILENAME
                        specify logfile name (default: sincho.log)
  -n NUMBER_CAND, --number_cand NUMBER_CAND
                        specify number of candidates. default: 10
```
**For the default execution, SINCHO is performed by the following command:**
```
$ sincho -p ${protein}.pdb -l ${ligand}.pdb
```
```${protein}.pdb``` and ```${ligand}.pdb``` must correspond to the input files in P2C.

## Output Files
output files in SINCHO are stored in current directory or ```./sincho-output/```.
In ```./sincho-output/```, ```cluster*.pqr.pdb``` is temporary file. ```check_terms.txt``` stored the Extend Score decomposition of each pair.
```./sincho.log``` stores the results of the ranking pairs based on the Extend Score.

## Visualization
The SINCHO results can be viewed if the process was successful by the following code:
```
$ pymol ./sincho-output/sincho.pse
```

