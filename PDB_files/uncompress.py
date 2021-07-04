import os
for item in os.listdir("."):
    if item.endswith("gz"):
        os.system(f"gunzip -f {item}")