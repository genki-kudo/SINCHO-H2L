#!/bin/sh

if [ -z "$3" ]
then
   echo :Usage ./$( basename $0 ) pdbid lig cpu
   exit
fi

#$1->inputのcomplex, $2->対象リガンドの残基名
comp=$1
lig=$2
cpu=$3
cwd=$(pwd)

#リガンドのパラメータ作成
grep ${lig} ${comp}.pdb > ${lig}.pdb
antechamber -i ${lig}.pdb -fi pdb -o ${lig}.prep -fo prepi -c bcc -at gaff -nc 0
parmchk2 -i ${lig}.prep -o ${lig}.frcmod -f prepi -s gaff

sed '/ H/d' ${comp}.pdb > ${comp}_noh.pdb
sed '/CONECT/d' ${comp}_noh.pdb > lala.pdb
mv lala.pdb ${comp}_noh.pdb

#複合体パラメータ作成
cat << EOF > tleap_in.txt
source leaprc.protein.ff14SB   #sourceで力場指定
source leaprc.gaff
source leaprc.water.tip3p
loadAmberPrep ./${lig}.prep   #loadPrep/paramsで.prep/.frcmod指定
loadamberparams ./${lig}.frcmod
protein = loadpdb ./${comp}_noh.pdb  #complex.pdbをproteinと命名   
addions2 protein Na+ 0  #addions2で総電価0に
addions2 protein Cl- 0
solvateBox protein TIP3PBOX 15.0  #箱づくりと水入れ
saveamberparm protein ./complex_wat.prmtop ./complex_wat.inpcrd  #.prmtopと.inpcrd生成
savepdb protein complex_wat.pdb  #複合体構造(準備済みの
check protein
quit
EOF

tleap -f tleap_in.txt

cat << EOF > ./convert.py
import parmed as pmd
amber = pmd.load_file('complex_wat.prmtop',xyz='complex_wat.inpcrd')
amber.save('complex_wat.gro')
amber.save('complex_wat.top')
quit()

EOF

#python3 convert.py
python convert.py

rm ANTECHAMBER*

#Energy-Minimization
cat << eof0 > ./min.mdp
; minim.mdp - used as input into grompp to generate em.tpr
integrator      = steep         ; Algorithm (steep = steepest descent minimization)
emtol           = 1000.0        ; Stop minimization when the maximum force < 100.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size
nsteps          = 10000         ; Maximum number of (minimization) steps to perform
; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist             = 1             ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type             = grid              ; Method to determine neighbor list (simple, grid)
coulombtype         = PME               ; Treatment of long range electrostatic interactions
rcoulomb            = 1.0               ; Short-range electrostatic cut-off
rvdw                = 1.0               ; Short-range Van der Waals cut-off
pbc                     = xyz           ; Periodic Boundary Conditions (yes/no)
eof0
rm \#*
gmx_mpi grompp -f ./min.mdp -c ./complex_wat.gro -p ./complex_wat.top -o ./em.tpr -maxwarn 1
gmx_mpi mdrun -v -deffnm ./em
# gmx grompp -f ./min.mdp -c ./complex_wat.gro -p ./complex_wat.top -o ./em.tpr -maxwarn 1
# gmx mdrun -v -deffnm ./em
rm \#*


#NVT
cat << eof1 > nvt.mdp
title           = complex_wat ff14SB NVT equilibration
define          = -DPOSRES      ; position restrain the protein
; Run parameters
integrator      = md            ; leap-frog integrator
nsteps          = 50000                ; 1 * 50000 = 500 ps
dt                  = 0.001            ; 1 fs
; Output control
nstxout         = 0             ; save coordinates every 1.0 ps
nstvout         = 0             ; save velocities every 1.0 ps
nstenergy       = 0             ; save energies every 1.0 ps
nstlog          = 5000          ; update log file every 1.0 ps
nstxout-compressed  = 5000       ; save compressed coordinates every 10.0 ps
; Bond parameters
continuation            = no            ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints                 = H-bonds   ; all bonds (even heavy atom-H bonds) constrained
lincs_iter                  = 1             ; accuracy of LINCS
lincs_order                 = 4             ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type             = grid              ; search neighboring grid cells
nstlist             = 10                ; 20 fs, largely irrelevant with Verlet
rcoulomb            = 1.0               ; short-range electrostatic cutoff (in nm)
rvdw                = 1.0               ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype         = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order           = 4         ; cubic interpolation
fourierspacing  = 0.16  ; grid spacing for FFT
; Temperature coupling is on
tcoupl          = V-rescale                 ; modified Berendsen thermostat
tc-grps         = Protein Non-Protein   ; two coupling groups - more accurate
tau_t           = 0.1     0.1           ; time constant, in ps
ref_t           = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl          = no            ; no pressure coupling in NVT
; Periodic boundary conditions
pbc             = xyz               ; 3-D PBC
; Dispersion correction
DispCorr        = EnerPres      ; account for cut-off vdW scheme
; Velocity generation
gen_vel         = no            ; assign velocities from Maxwell distribution
gen_temp        = 100           ; temperature for Maxwell distribution
gen_seed        = -1            ; generate a random seed
eof1
#comm-mode=LINEAR
#comm-grps=Protein
#nstcomm=1
#nstcalcenergy=1
rm \#*
gmx_mpi grompp -maxwarn 3 -f nvt.mdp -c em.gro -p complex_wat.top -o nvt.tpr -r em.gro
# gmx grompp -maxwarn 3 -f nvt.mdp -c em.gro -p complex_wat.top -o nvt.tpr -r em.gro
cwd=$(pwd)
cd ${cwd}
OMP_NUM_THREADS=${cpu}
gmx_mpi mdrun -deffnm nvt -ntomp ${cpu}
# gmx mdrun -deffnm nvt #-ntomp ${cpu} 
eof2
qsub qsub_nvt.sh
rm \#*

#NPT
cat << eof1 > ./npt.mdp
title           = complex_wat NPT equilibration
define          = -DPOSRES      ; position restrain the protein
; Run parameters
integrator      = md            ; leap-frog integrator
nsteps          = 50000         ; 1 * 50000 = 500 ps
dt                  = 0.001             ; 1 fs
; Output control
nstxout         = 0             ; save coordinates every 1.0 ps
nstvout         = 0             ; save velocities every 1.0 ps
nstenergy       = 0             ; save energies every 1.0 ps
nstlog          = 500           ; update log file every 1.0 ps
nstxout-compressed  = 500       ; save compressed coordinates every 10.0 ps
; Bond parameters
continuation            = yes           ; Restarting after NVT
constraint_algorithm    = lincs     ; holonomic constraints
constraints                 = all-bonds ; all bonds (even heavy atom-H bonds) constrained
lincs_iter                  = 1             ; accuracy of LINCS
lincs_order                 = 4             ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type             = grid              ; search neighboring grid cells
nstlist             = 10            ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb            = 1.0               ; short-range electrostatic cutoff (in nm)
rvdw                = 1.0               ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype         = PME               ; Particle Mesh Ewald for long-range electrostatics
pme_order           = 4             ; cubic interpolation
fourierspacing  = 0.16          ; grid spacing for FFT
; Temperature coupling is on
tcoupl          = V-rescale                 ; modified Berendsen thermostat
tc-grps         = Protein Non-Protein   ; two coupling groups - more accurate
tau_t           = 0.1     0.1           ; time constant, in ps
ref_t           = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Berendsen         ; Pressure coupling on in NPT
pcoupltype              = isotropic                 ; uniform scaling of box vectors
tau_p                   = 2.0                       ; time constant, in ps
ref_p                   = 1.0                       ; reference pressure, in bar
compressibility     = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com
; Periodic boundary conditions
pbc             = xyz           ; 3-D PBC
; Dispersion correction
DispCorr        = EnerPres      ; account for cut-off vdW scheme
; Velocity generation
gen_vel         = no            ; Velocity generation is off
eof1
#comm-mode=LINEAR
#comm-grps=Protein
#nstcomm=1
#nstcalcenergy=1
rm \#*
gmx_mpi grompp -maxwarn 3 -f npt.mdp -c nvt.gro -p complex_wat.top -o npt.tpr -r nvt.gro
# gmx grompp -maxwarn 3 -f npt.mdp -c nvt.gro -p complex_wat.top -o npt.tpr -r nvt.gro
cwd=$(pwd)
cd ${cwd}
OMP_NUM_THREADS=${cpu}
gmx_mpi mdrun -deffnm npt -ntomp ${cpu}
# gmx mdrun -deffnm npt #-ntomp ${cpu} 
rm \#*

cp npt.gro prod0.gro

#Production-Run
id=1
pre=`expr $id - 1`
cat << eof1 > prod.mdp
title           = complex_wat MD simulation (NPT)
;define          = -DPOSRES      ; position restrain the protein
; Run parameters
integrator      = md            ; leap-frog integrator
nsteps          = 5000              ; 2 * 5000 = 10 ps
dt                  = 0.002             ; 2 fs
; Output control
nstxout         = 0           ; save coordinates every 1.0 ps
nstvout         = 0           ; save velocities every 1.0 ps
nstenergy       = 0           ; save energies every 1.0 ps
nstlog          = 500           ; update log file every 50 ps
nstxout-compressed = 500
; Bond parameters
continuation            = no           ; Restarting after NVT
constraint_algorithm    = lincs     ; holonomic constraints
constraints                 = all-bonds ; all bonds (even heavy atom-H bonds) constrained
lincs_iter                  = 1             ; accuracy of LINCS
lincs_order                 = 4             ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type             = grid              ; search neighboring grid cells
nstlist             = 10            ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb            = 1.0             ; short-range electrostatic cutoff (in nm)
rvdw                = 1.0               ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype         = PME               ; Particle Mesh Ewald for long-range electrostatics
pme_order           = 4             ; cubic interpolation
fourierspacing  = 0.16          ; grid spacing for FFT
; Temperature coupling is on
tcoupl          = V-rescale                 ; modified Berendsen thermostat
tc-grps         = Protein Non-Protein   ; two coupling groups - more accurate
tau_t           = 0.1     0.1           ; time constant, in ps
ref_t           = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Berendsen         ; Pressure coupling on in NPT
pcoupltype              = isotropic                 ; uniform scaling of box vectors
tau_p                   = 2.0                       ; time constant, in ps
ref_p                   = 1.0                       ; reference pressure, in bar
compressibility     = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com
; Periodic boundary conditions
pbc             = xyz           ; 3-D PBC
; Dispersion correction
DispCorr        = EnerPres      ; account for cut-off vdW scheme
; Velocity generation
gen_vel         = yes            ; Velocity generation is off
gen_seed  = -1
eof1
gmx_mpi grompp -maxwarn 3 -f prod.mdp -c prod$pre.gro -r prod$pre.gro -p complex_wat.top -o prod$id.tpr
# gmx grompp -maxwarn 3 -f prod.mdp -c prod$pre.gro -r prod$pre.gro -p complex_wat.top -o prod$id.tpr
OMP_NUM_THREADS=${cpu}
gmx_mpi mdrun -deffnm prod${id} -ntomp ${cpu}
# gmx mdrun -deffnm prod${id} #-ntomp ${cpu} 

rm \#*
# gmx make_ndx -f complex_wat.gro <<EOF
gmx_mpi make_ndx -f complex_wat.gro <<EOF
del16
del16
!19
name 20 nowation
q
EOF

# gmx trjconv -f prod1.xtc -o prod1nj.xtc -s em.tpr -n index.ndx -pbc whole <<EOF
gmx_mpi trjconv -f prod1.xtc -o prod1nj.xtc -s em.tpr -n index.ndx -pbc whole <<EOF
System
EOF

# gmx trjconv -f prod1nj.xtc -o prod1nj2.xtc -s em.tpr -n index.ndx -pbc cluster <<EOF
gmx_mpi trjconv -f prod1nj.xtc -o prod1nj2.xtc -s em.tpr -n index.ndx -pbc cluster <<EOF
nowation
System
EOF

# gmx trjconv -f prod1nj2.xtc -o prod1nj3.xtc -s em.tpr -n index.ndx -pbc mol -ur compact -center <<EOF
gmx_mpi trjconv -f prod1nj2.xtc -o prod1nj3.xtc -s em.tpr -n index.ndx -pbc mol -ur compact -center <<EOF
nowation
System
EOF

# gmx trjconv -f prod1nj3.xtc -o prod1nj4.xtc -s em.tpr -n index.ndx -center -fit rot+trans <<EOF
gmx_mpi trjconv -f prod1nj3.xtc -o prod1nj4.xtc -s em.tpr -n index.ndx -center -fit rot+trans <<EOF
nowation
nowation
System
EOF

# gmx trjconv -f prod1nj4.xtc -s complex_wat.pdb -o conformation1-10.pdb -center -pbc whole -n index.ndx <<EOF
gmx_mpi trjconv -f prod1nj4.xtc -s complex_wat.pdb -o conformation1-10.pdb -center -pbc whole -n index.ndx <<EOF
C-alpha
nowation
EOF





