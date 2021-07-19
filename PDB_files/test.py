from pymol import cmd
cmd.load("2wmb.pdb")
cmd.select("/2wmb//A")
cmd.extract("A", "sele")
cmd.deselect()
cmd.select("/2wmb//C")
cmd.extract("C", "sele")
cmd.deselect()
cmd.align("C","A")
cmd.delete("2wmb")
cmd.deselect()
cmd.remove("resn HOH")

cmd.alter("all", "segi=''")
cmd.alter("all", "chain=''")
print(f"altered")
cmd.deselect()
cmd.select("/A///0:4/CA")
#/A///0:4/CA
#/B///0:4/CA
#rms_cur /A///0:4/CA, /B///0:4/CA

xyz = cmd.get_coords('sele', 1) 
print(f"xyz: {xyz}")

cmd.deselect()
cmd.select("A///199/CA")
xyz199 = cmd.get_coords('sele', 1) 
print(f"xyz199: {xyz199}")