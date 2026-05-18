import os
import shutil

ff53_dir = r"c:\Users\cores\NoOneDrive\Projects\Increment-mitomycin_c-bioavailability\DPPC\gromos53a6_lipid.ff"
ff54_dir = r"c:\Users\cores\NoOneDrive\Projects\Increment-mitomycin_c-bioavailability\mitomycin_c\gromos54a7_atb.ff"

ff53_nb = os.path.join(ff53_dir, "ffnonbonded.itp")
ff54_nb = os.path.join(ff54_dir, "ffnonbonded.itp")

# Create a backup
if not os.path.exists(ff53_nb + ".bak"):
    shutil.copy(ff53_nb, ff53_nb + ".bak")

def parse_itp(file_path):
    data = {"[ atomtypes ]": [], "[ nonbond_params ]": [], "[ pairtypes ]": []}
    if not os.path.exists(file_path): return data
    with open(file_path, "r") as f:
        lines = f.readlines()
    current_sec = None
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("["):
            sec = stripped.split("]")[0] + "]"
            if sec in data:
                current_sec = sec
            else:
                current_sec = None
            continue
        if current_sec and not stripped.startswith(";") and len(stripped) > 0:
            data[current_sec].append(line)
    return data

data53 = parse_itp(ff53_nb)
data54 = parse_itp(ff54_nb)

# Find existing
def get_sig(line, sec):
    parts = line.split()
    if sec == "[ atomtypes ]":
        return parts[0] if len(parts)>0 else None
    elif sec == "[ nonbond_params ]":
        if len(parts) >= 2:
            return "-".join(sorted([parts[0], parts[1]]))
    elif sec == "[ pairtypes ]":
        if len(parts) >= 2:
            return "-".join(sorted([parts[0], parts[1]]))
    return None

existing_keys = {
    "[ atomtypes ]": set(get_sig(l, "[ atomtypes ]") for l in data53["[ atomtypes ]"] if get_sig(l, "[ atomtypes ]")),
    "[ nonbond_params ]": set(get_sig(l, "[ nonbond_params ]") for l in data53["[ nonbond_params ]"] if get_sig(l, "[ nonbond_params ]")),
    "[ pairtypes ]": set(get_sig(l, "[ pairtypes ]") for l in data53["[ pairtypes ]"] if get_sig(l, "[ pairtypes ]"))
}

to_append = {
    "[ atomtypes ]": [],
    "[ nonbond_params ]": [],
    "[ pairtypes ]": []
}

for sec in ["[ atomtypes ]", "[ nonbond_params ]", "[ pairtypes ]"]:
    for line in data54[sec]:
        sig = get_sig(line, sec)
        if sig and sig not in existing_keys[sec]:
            to_append[sec].append(line)

print(f"Missing atomtypes: {len(to_append['[ atomtypes ]'])}")
print(f"Missing nonbond_params: {len(to_append['[ nonbond_params ]'])}")
print(f"Missing pairtypes: {len(to_append['[ pairtypes ]'])}")

# Now read the original file line by line and append right before the next section
with open(ff53_nb + ".bak", "r") as f:
    lines53 = f.readlines()

new_lines = []
current_sec = None
appended = {"[ atomtypes ]": False, "[ nonbond_params ]": False, "[ pairtypes ]": False}

for line in lines53:
    stripped = line.strip()
    if stripped.startswith("["):
        sec = stripped.split("]")[0] + "]"
        if current_sec in appended and not appended[current_sec]:
            new_lines.append(f"\n; --- Appended from ATB 54a7 {current_sec} ---\n")
            new_lines.extend(to_append[current_sec])
            new_lines.append("\n")
            appended[current_sec] = True
        current_sec = sec
    new_lines.append(line)

# Handle the last section
if current_sec in appended and not appended[current_sec]:
    new_lines.append(f"\n; --- Appended from ATB 54a7 {current_sec} ---\n")
    new_lines.extend(to_append[current_sec])
    new_lines.append("\n")
    appended[current_sec] = True

with open(ff53_nb, "w") as f:
    f.writelines(new_lines)
print("ffnonbonded.itp updated.")

# 2. Update atomtypes.atp
ff53_atp = os.path.join(ff53_dir, "atomtypes.atp")
ff54_atp = os.path.join(ff54_dir, "atomtypes.atp")
if os.path.exists(ff53_atp) and os.path.exists(ff54_atp):
    if not os.path.exists(ff53_atp + ".bak"):
        shutil.copy(ff53_atp, ff53_atp + ".bak")
    with open(ff53_atp + ".bak", "r") as f:
        lines53_atp = f.readlines()
    atp_types53 = set(line.split()[0] for line in lines53_atp if len(line.split()) > 1 and not line.startswith(";"))
    
    with open(ff54_atp, "r") as f:
        lines54_atp = f.readlines()
    atp_append = [line for line in lines54_atp if len(line.split()) > 1 and not line.startswith(";") and line.split()[0] not in atp_types53]
    
    with open(ff53_atp, "w") as f:
        f.writelines(lines53_atp)
        f.write("\n; --- Appended ATB atomtypes from 54a7 ---\n")
        f.writelines(atp_append)
    print(f"Added {len(atp_append)} atom types to atomtypes.atp")

print("Done blending 54A7 ATB parameters into 53A6 lipid force field.")
