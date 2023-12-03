# SINCHO(Site Idenrification and Next CHOice) protocol for Hit-to-Lead

SINCHO protocol is the method for **the prediction/suggestion of the desirable anchor atom and growth site pair** for the modification of the hitcompound in hit-to-lead process.

##Installation
###Requirements
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
  In the SINCHO protocol, you should additionally install below modules.
  * Biopython(1.79)
    ```
    conda install Biopython
  


  
##Input files
* Ligand pdb file
  *水素付加済
  *原子名付
  *ファイル名=残基名
  *CONECT行付
