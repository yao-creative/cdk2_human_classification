from pymol import cmd
cmd.load("1dm2.pdb")
cmd.select("/1dm2//A")
cmd.extract("A", "sele")
"""cmd.deselect()
cmd.select("/1dm2//C")
cmd.extract("C", "sele")
cmd.deselect()
cmd.align("C","A")
cmd.delete("1dm2")"""
cmd.deselect()
cmd.load("2wih.pdb")
cmd.select("/2wih//A")
cmd.extract("A2", "sele")
cmd.deselect()
cmd.select("/2wih//C")
cmd.extract("C2", "sele")
cmd.deselect()
cmd.align("C2","A")
cmd.align("A2","A")
cmd.delete("2wih")
cmd.deselect()
cmd.remove("resn HOH")
list = cmd.id_atom("/A///2/CA")
print(f"id_atom: {list}")
cmd.alter("all", "segi=''")
cmd.alter("all", "chain=''")
print(f"altered")
cmd.deselect()
cmd.select("/A///0:4/CA")
#/A///0:4/CA
#/B///0:4/CA
# 1DM2_A other_id: 2WIH_A
#rms_cur /A///0:4/CA, /B///0:4/CA

xyz = cmd.get_coords('sele', 1) 
print(f"xyz: {xyz}")

cmd.deselect()
"""cmd.select("/A///1:6/CA")
cmd.create("A_aln", "sele")
cmd.deselect()
cmd.select("/C///1:6/CA")
cmd.create("C_aln", "sele")
cmd.deselect()"""
rms= cmd.rms_cur("/A////CA","/A2////CA",matchmaker=4)

#rms=cmd.rms_cur("A_aln","C_aln",matchmaker=-1)
print(f"rms found! {rms}")