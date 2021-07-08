import os
#bash calls a python script which immitates bash : D
for item in os.listdir("."):
    if item.endswith("gz"):
        os.system(f"gunzip -f {item}")