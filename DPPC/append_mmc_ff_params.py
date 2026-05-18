import os
import shutil

ff53_dir = r"c:\Users\cores\NoOneDrive\Projects\Increment-mitomycin_c-bioavailability\DPPC\gromos53a6_lipid.ff"
ff54_dir = r"c:\Users\cores\NoOneDrive\Projects\Increment-mitomycin_c-bioavailability\mitomycin_c\gromos54a7_atb.ff"

# 1. Update ffnonbonded.itp
ff53_nb = os.path.join(ff53_dir, "ffnonbonded.itp")
ff54_nb = os.path.join(ff54_dir, "ffnonbonded.itp")

# Create a backup
shutil.copy(ff53_nb, ff53_nb + ".bak")

# Read 53a6 types
with open(ff53_nb, "r") as f:
    lines53_nb = f.readlines()
types53 = set()
in_atomtypes = False
for line in lines53_nb:
    if line.startswith("[ atomtypes ]"):
        in_atomtypes = True
    elif line.startswith("[") and not line.startswith("[ atomtypes ]"):
        in_atomtypes = False
    if in_atomtypes and len(line.split()) > 5 and not line.startswith(";"):
        types53.add(line.split()[0])

# Read 54a7 types and append missing
with open(ff54_nb, "r") as f:
    lines54_nb = f.readlines()
in_atomtypes = False
to_append = []
for line in lines54_nb:
    if line.startswith("[ atomtypes ]"):
        in_atomtypes = True
        continue
    elif line.startswith("["):
        in_atomtypes = False
        continue
    
    if in_atomtypes and len(line.split()) > 5 and not line.startswith(";"):
        typ = line.split()[0]
        if typ not in types53:
            to_append.append(line)

# Insert the missing types right before the next directive or at the end of [ atomtypes ]
new_lines53_nb = []
inserted = False
for line in lines53_nb:
    if line.startswith("[ nonbond_params ]") and not inserted:
        new_lines53_nb.append("; --- Appended ATB atomtypes from 54a7 ---\n")
        new_lines53_nb.extend(to_append)
        new_lines53_nb.append("\n")
        inserted = True
    new_lines53_nb.append(line)

if not inserted: # fallback
    new_lines53_nb.append("\n; --- Appended ATB atomtypes from 54a7 ---\n")
    new_lines53_nb.extend(to_append)    

with open(ff53_nb, "w") as f:
    f.writelines(new_lines53_nb)
print(f"Added {len(to_append)} atom types to ffnonbonded.itp")


# 2. Update atomtypes.atp (if it exists)
ff53_atp = os.path.join(ff53_dir, "atomtypes.atp")
ff54_atp = os.path.join(ff54_dir, "atomtypes.atp")
if os.path.exists(ff53_atp) and os.path.exists(ff54_atp):
    shutil.copy(ff53_atp, ff53_atp + ".bak")
    with open(ff53_atp, "r") as f:
        lines53_atp = f.readlines()
    atp_types53 = set(line.split()[0] for line in lines53_atp if len(line.split()) > 1 and not line.startswith(";"))
    
    with open(ff54_atp, "r") as f:
        lines54_atp = f.readlines()
    atp_append = [line for line in lines54_atp if len(line.split()) > 1 and not line.startswith(";") and line.split()[0] not in atp_types53]
    
    with open(ff53_atp, "a") as f:
        f.write("\n; --- Appended ATB atomtypes from 54a7 ---\n")
        f.writelines(atp_append)
    print(f"Added {len(atp_append)} atom types to atomtypes.atp")

print("Done blending 54A7 ATB parameters into 53A6 lipid force field.")
