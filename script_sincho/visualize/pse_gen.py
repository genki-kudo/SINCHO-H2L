from pymol import cmd
import glob
import os

def visualize(protein, ligand, check_terms_sorted, poc_clus_dir, outdir, num):
    cmd.load(ligand)
    cmd.load(protein)
    cmd.color("atomic", ligand.split("/")[-1][:-4])
    cmd.color("gray70", protein.split("/")[-1][:-4])
    clsname = glob.glob(poc_clus_dir+'/*.pqr')
    for n in range(1,len(clsname)+1):
        cls = poc_clus_dir+'/cluster'+str(n)+'.pqr'
        cmd.load(cls)
        number = os.path.splitext(os.path.basename(cls))[0][7:]
        cmd.color(int(number)+2, os.path.splitext(os.path.basename(cls))[0]) 
        cmd.show("mesh", os.path.splitext(os.path.basename(cls))[0])
        cmd.show("spheres", os.path.splitext(os.path.basename(cls))[0])
        cmd.hide("sticks", os.path.splitext(os.path.basename(cls))[0])
    cmd.set("sphere_scale", "0.3")
    cmd.util.cbag(ligand.split("/")[-1][:-4])
    cmd.bg_color("white")
    cmd.zoom("all")
    cmd.scene("p2c_result", "store")
    cmd.save(outdir+"/sincho.pse")

    for i in open(ligand):
        if i[0:6]=="ATOM  "or i[0:6]=="HETATM":
            ligres = i[17:20]
            break

    r=0
    for n in check_terms_sorted:
    #for i in open(check_terms_sorted):
        #n = i.split()

        if r<=int(num):
            rank = "rank"+str(r+1)
            cls = n[0][:-4]
            #atom = ligand.split("/")[-1][:-4]+"/"+n[1].split("_")[1]
            atom = ligres+"/"+n[1].split("_")[1]
            print(atom)
            score = n[5]
            cmd.color("gray70", "all")
            cmd.select(atom)
            cmd.create(rank+"_atom", "sele")
            #cmd.color("atomic", ligand.split("/")[-1][:-4])
            cmd.util.cbag(ligand.split("/")[-1][:-4])
            cmd.color("red", cls)
            cmd.color("red", rank+"_atom")
            cmd.show("spheres", rank+"_atom")
            cmd.distance(rank+"_mark", rank+"_atom", cls, mode=4)
            cmd.hide("label",rank+"_mark")
            cmd.color("red",rank+"_mark")
            cmd.set("dash_width","5.0")
            cmd.group(rank, rank+"_atom "+rank+"_mark")
            cmd.select(cls+"|"+ligand.split("/")[-1][:-4])
            cmd.orient("sele")
            cmd.scene(rank+"_"+str(score), "store")
            cmd.hide("spheres",rank+"_atom")
            cmd.hide("dash",rank+"_mark")
        r+=1
    
    cmd.zoom("all")
    cmd.save(outdir+"/sincho.pse")
    cmd.delete("all")


