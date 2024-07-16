from decimal import Decimal
import parmed as pmd 
import copy

def read_top_file(top_file_name):
    top_file = open(top_file_name, "r")
    allowed_keywords = ["atoms", "atomtypes", "bonds", "angles", "dihedrals"]
    data = []
    current_data_dict = {}
    current_keyword = ""
    for line in top_file:
        if line.startswith("["):
            keyword = line.strip("[ ").strip(" ]\n")
            current_keyword = keyword
            if keyword == "moleculetype":
                data.append(copy.deepcopy(current_data_dict))
                current_data_dict = {}
            if keyword == "system":
                data.append(copy.deepcopy(current_data_dict))
                break
            continue
        if line.startswith(";") or line == "\n" or line.strip(" ") == "\n" or line.startswith("#") or current_keyword not in allowed_keywords:
            continue 
        else:
            if current_keyword not in current_data_dict.keys():
                    current_data_dict[current_keyword] = [getVals(line,";")]
            else:
                    current_data_dict[current_keyword].append(getVals(line, ";"))
    atomtype_list = []
    
    for i in range(0,len(data)):
        date = data[i]
        current_type_dict = {}
        for key in date.keys():
            dat_full = date[key]
            if key == "atomtypes":
                for dat in dat_full:
                    atomtype_list.append([dat[0],dat[2],dat[5],dat[6]])
                continue 
            if key == "atoms":
                for dat in dat_full:
                    current_type_dict[dat[0]] = dat[1]
                data[i][key] = current_type_dict
            if key == "bonds":  
                bonds_fixed = []
                for bond in dat_full:
                    bond[0] = current_type_dict[bond[0]]
                    bond[1] = current_type_dict[bond[1]]
                    bonds_fixed.append(bond)
                data[i][key] = bonds_fixed
            if key == "angles":
                angles_fixed = []
                for angle in dat_full:
                    angle[0] = current_type_dict[angle[0]]
                    angle[1] = current_type_dict[angle[1]]
                    angle[2] = current_type_dict[angle[2]]
                    angles_fixed.append(angle)
                data[i][key] = angles_fixed
            if key == "dihedrals":
                dihedrals_fixed = []
                for dihedral in dat_full:
                    dihedral[0] = current_type_dict[dihedral[0]]
                    dihedral[1] = current_type_dict[dihedral[1]]
                    dihedral[2] = current_type_dict[dihedral[2]]
                    dihedral[3] = current_type_dict[dihedral[3]]
                    dihedrals_fixed.append(dihedral)
                data[i][key] = dihedrals_fixed
    data.pop(0)
    top_file.close()
    return atomtype_list, data

def replace_angle(line_in, angle_pre, angle_now):
    angle_len = len(angle_pre)
    angle_pos = line_in.find(angle_pre)
    res = line_in[:angle_pos] + angle_now + line_in[angle_pos+angle_len:]
    return res

def getVals(line_in, end_sym):
    line_sim = line_in.find(end_sym)
    if line_sim != -1:
        line_in = line_in[0:line_sim]
    line_in = line_in.strip("\n").split(" ")
    data = []
    for a in line_in:
        if a != "" and a != "\t" and a != " ":
            data.append(a) 
    return data

def get_permutations_3(in_list):
    perm = [in_list[0]+in_list[1]+in_list[2],in_list[2]+in_list[1]+in_list[0]]
    for i in range(0,3):
        for j in range(0,3):
            for k in range(0,3):
                if i!= j and j != k and k != i:
                    perm.append(in_list[i]+in_list[j]+in_list[k])
    
    return setO(perm)

def get_permutations_4(in_list):
    perm = [in_list[0]+in_list[1]+in_list[2]+in_list[3], in_list[3]+in_list[2]+in_list[1]+in_list[0]]
    for i in range(0,4):
        for j in range(0,4):
            for k in range(0,4):
                for m in range(0,4):
                    if i != j and i != k and i != m and j != k and j!= m and k!= m:
                        perm.append(in_list[i]+in_list[j]+in_list[k]+in_list[m])
    return setO(perm)

def get_permutations_core(in_list):
    perm = [in_list[0]+in_list[1]+in_list[2]+in_list[3], in_list[3]+in_list[2]+in_list[1]+in_list[0]]
    perm.append(in_list[0]+in_list[1]+in_list[2])
    perm.append(in_list[1]+in_list[2]+in_list[3])
    perm.append(in_list[2]+in_list[1]+in_list[0])
    perm.append(in_list[3]+in_list[2]+in_list[1])
    perm.append(in_list[1]+in_list[2])
    perm.append(in_list[2]+in_list[1])
    return perm 
def get_permutations_4_wild(in_list):
    perm = []
    perm = list(get_permutations_4(in_list))
    temp_a = perm[2]
    temp_b = perm[3]
    perm[2] = in_list[1]+in_list[2]
    perm[3] = in_list[2]+in_list[1]
    perm.append(temp_a)
    perm.append(temp_b)
    for i in range(0,4):
        for j in range(0,4):
            for k in range(0,4):
                if i != j and i != k and j != k:
                    perm.append(in_list[i]+in_list[j]+in_list[k])
    for i in range(0,4):
        for j in range(0,4):
            if i!=j:
                perm.append(in_list[i]+in_list[j])
    return setO(perm)

def parse_inputline(line_in):
    dat_pre = parse_pre(getVals(line_in, "/"))
    dat_pre = dat_pre[0:4]  
    bond_type = dat_pre[0] 
    data = [dat_pre[0], dat_pre[1].strip("-"), bondDistConv(dat_pre[3].strip("\t")), forceFactorConv(dat_pre[2].strip("\t"))]
    return data

def parse_inputline_angle(line_in):
    data_pre = parse_pre(getVals(line_in, "/"))
    temp = str(round(Decimal(data_pre[3])*Decimal(2)*Decimal(4.184),6))
    data_pre[3] = data_pre[4]
    data_pre[4] = temp
    return data_pre

def parse_inputline_dihedral_prop(line_in):
    data_pre = parse_pre(getVals(line_in, "SCEE"))  
    if data_pre == []:
        return data_pre
    temp = str(round(float(data_pre[5])*float(4.184)/float(data_pre[4]),7))
    data_pre = data_pre[:8]
    data_pre[4] = 1
    data_pre[5] = data_pre[6]
    data_pre[6] = temp
    data_pre[7] = str(abs(float(data_pre[7])))
    return data_pre

def parse_inputline_dihedral_improp(line_in):
    data_pre = parse_pre(getVals(line_in, "SCEE"))
    data_res = [0,0,0,0,0,0,0,0]
    for i in range(0,4):
        data_res[i] = data_pre[i]
    data_res[5] = data_pre[5]
    data_res[6] = str(round(float(data_pre[4])*float(4.184),7))
    data_res[7] = data_pre[6]
    data_res[4] = 4
    return data_res


def findMatch(dihedral, dihedral_list):
    for i in range(0,len(dihedral_list)):
        if float(dihedral[6]) == float(dihedral_list[i][2]):
            return i
    return -1

def parse_pre(data_preelim):
    res = []
    for a in data_preelim:
        a = a.strip(" ").strip("\t")
        try:
            float(a)
            ar = [a]
        except:
            a = a.strip("-")
            ar = a.split("-") 
        
        for aar in ar: 
            res.append(aar) 
    return res


def convertR(inV):
    return str(Decimal(inV)/(Decimal(5)*Decimal(2)**(Decimal(1)/Decimal(6))))
def convertEps(inEps):
    return str(round(Decimal(inEps)*Decimal(4.184),8))

def get_file_no_suffix(file_name):
    pos = -1
    while True:
        if file_name.find(".", pos+1) == -1:
            break
        else:
            pos = file_name.find(".",pos+1)
    if pos == -1:
        print("ERROR: No file suffix in secondary input file")
    return file_name[0:pos]

def fix_secondary_file(file_in_name):
    sec_input_file = open(file_in_name, "r") 
    file_in_name_fixed = get_file_no_suffix(file_in_name) + "_fixed.dat"
    sec_input_file_fixed = open(file_in_name_fixed, "w")
    first_mass = True
    first_angle = True
    first_bond = True
    first_dihedral_prop = True
    first_dihedral_improp = True

    for line in sec_input_file:
        current_line = line_det_type(line)
        if current_line == 1 and first_mass:
            first_mass = False
            sec_input_file_fixed.write("MASS\n")
        if current_line == 2 and first_bond:
            first_bond = False
            sec_input_file_fixed.write("BOND\n")
        if current_line == 3 and first_angle:
            first_angle = False
            sec_input_file_fixed.write("ANGLE\n")
        if current_line == 4 and first_dihedral_prop:
            first_dihedral_prop = False
            sec_input_file_fixed.write("DIHEDRAL\n")
        if current_line == 5 and first_dihedral_improp:
            first_dihedral_improp = False
            sec_input_file_fixed.write("IMPROPER\n")
        sec_input_file_fixed.write(line)

    sec_input_file.close()
    sec_input_file_fixed.close()
    return file_in_name_fixed

def line_det_type(line):
    dat = parse_pre(getVals(line,";"))
    first_atom = False
    second_atom = False
    third_atom = False
    fourth_atom = False
    # 0 Empty line, 1 Mass, 2 Bond, 3 Angle, 4 Dihedral prop, 5 Dihedral improp, 6 Something else
    if len(dat) < 3: 
        return 0
    try:
        float(dat[0])
    except:
        first_atom = True
    try: 
        float(dat[1])
    except: 
        second_atom = True

    try:
        float(dat[2])
    except: 
        third_atom = True

    try:
        float(dat[3])
    except: 
        fourth_atom = True

    if first_atom and second_atom and third_atom and fourth_atom:
        try: 
            float(dat[4])
            try:
                float(dat[5])
                float(dat[6])
                float(dat[7])
                return 4
            except:
                return 5
        except:
            return 6
    if first_atom and second_atom and third_atom:
        return 3
    if first_atom and second_atom:
        return 2
    if first_atom:
        return 1

    return -1


def bondDistConv(bondDist):
    return str(Decimal(bondDist)/Decimal(10))

def forceFactorConv(forceFact):
    return str(round(Decimal(forceFact)*Decimal(100)*Decimal(2)*Decimal(4.184),4))

def setO(list_in):
    seen = set()
    return [x for x in list_in if not (x in seen or seen.add(x))]


one_step_prep = input("Do you want to create a .gro and .top file from a .prmtop and a .inpcrd file? (The generated files will be used in the angle correction process)[y/n]\n")
if one_step_prep == "y":
    prmtop_file_name = input("Provide the filename of the .prmtop file (with suffix)\n")
    inpcrd_file_name = input("Provide the filename of the .inpcrd file (with suffix)\n")
    input_file1_name = input("Provide the filename of the primary FF input file. (in case of tumuc this is tumuc.frcmod)\n")
    input_file2_name = input("Provide the filename of the secondary FF input file. (in case of tumuc this is parm10.dat)\n")
    top_base = get_file_no_suffix(prmtop_file_name)
    coord_base = get_file_no_suffix(inpcrd_file_name)
    parm = pmd.load_file(prmtop_file_name, inpcrd_file_name)
    parm.save(top_base+".top", format='gromacs')
    parm.save(coord_base + ".gro")
    top_file_name = top_base + ".top"
    gro_file_name = coord_base + ".gro" 
    print(f"Coordinates successfully written to {gro_file_name}.")
    print(f"Topology succesfully written to {top_file_name}.")
    print("\n")
else:
    input_wished = input("Do you want to manually select files for the forcefield check and possible angle correction? [y/n])\n")
    if input_wished == "y":
        input_file1_name = input("Provide the filename of the primary FF input file. (in case of tumuc this is tumuc.frcmod)\n")
        input_file2_name = input("Provide the filename of the secondary FF input file. (in case of tumuc this is parm10.dat)\n")
        top_file_name = input("Provide the filename of the GROMACS topology file.\n")
    else:
        # Defaults
        input_file1_name = "tumuc.frcmod"
        input_file2_name = "parm10.dat"
        top_file_name = "Abulge_fixed.top"


atomtype_list, data_raw = read_top_file(top_file_name)

input_file2_name = fix_secondary_file(input_file2_name)


print("\n")
print("Top file processed")
input_file = open(input_file1_name, "r") 

bonds = False 
angles = False
dihedrals_prop = False
dihedrals_improp = False
mass = False
nonbond = False
bond_input_dict = {}
angle_input_dict = {}
dihedral_prop_input_dict = {}
dihedral_improp_input_dict = {}
mass_nonbond_input_dict = {}
for line in input_file:
    if line.startswith("MASS"):
        mass = True
        continue 
    if line.startswith("BOND"):
        bonds = True
        continue
    if line.startswith("ANGLE"):
        angles = True
        continue 
    if line.startswith("DIHEDRAL"):
        dihedrals_prop = True
        continue
    if line.startswith("IMPROPER"):
        dihedrals_improp = True
        continue
    if line.startswith("NONBOND") or line.startswith("MOD4"):
        nonbond = True
        continue 
    if line == "\n" or line == " " or line.strip(" ") == "\n" or line == "\t":
        angles = False
        bonds = False
        dihedrals_prop = False
        dihedrals_improp = False
        mass = False
        nonbond = False
        continue
    if bonds:
        dat = parse_inputline(line)
        bond_input_dict[dat[0]+dat[1]] = dat[2:]
    if angles: 
        dat = parse_inputline_angle(line)
        angle_input_dict[dat[0]+dat[1]+dat[2]] = dat[3:]
    if dihedrals_prop:
        dat = parse_inputline_dihedral_prop(line)
        names = dat[0:4]
        name = ""
        for a in names:
            if a != "X":
                name+=a
        if name not in dihedral_prop_input_dict.keys():
            dihedral_prop_input_dict[name] = [dat[4:]]
        else:
            dihedral_prop_input_dict[name].append(dat[4:])
    if dihedrals_improp:
        dat = parse_inputline_dihedral_improp(line)
        names = dat[0:4]
        name = ""
        for a in names:
            if a != "X":
                name+=a
        if name not in dihedral_improp_input_dict.keys():
            dihedral_improp_input_dict[name] = [dat[4:]]
        else:
            dihedral_improp_input_dict[name].append(dat[4:])
    if mass:
        dat = getVals(line,"/")[:2]
        new_dat = []
        for a in dat:
            a = a.split("\t")
            for aa in a:
                if aa!="" and aa != " " and aa != "\n":
                    new_dat.append(aa)
        dat = new_dat
        mass_nonbond_input_dict[dat[0]] = [dat[1]]
    if nonbond: 
        dat = getVals(line, "/")[:3]
        if dat[0] in mass_nonbond_input_dict.keys():
            mass_nonbond_input_dict[dat[0]].append(convertR(dat[1]))
            mass_nonbond_input_dict[dat[0]].append(convertEps(dat[2]))
        else:
            print("Atom from nonbond parameters not found in mass directive")

print("Primary FF imput processed") 
input_file_2 = open(input_file2_name,"r")

bonds = False
dihedrals_prop = False
dihedrals_improp = False
angles = False 
mass = False
nonbond = False
cache = []
mass_nonbond_input_dict_parm = {}
dihedral_prop_input_dict_parm = {}
dihedral_improp_input_dict_parm = {}
bond_input_dict_parm = {}
angle_input_dict_parm = {} 
alias_dict = {}
for line in input_file_2:
    cache.append(line)
    if len(cache) > 10:
        cache = cache[len(cache)-10:len(cache)]
    if line.startswith("ANGLE"):
        angles = True
        continue 
    if line.startswith("DIHEDRAL"):
        dihedrals_prop = True
        continue
    if line.startswith("IMPROPER"):
        dihedrals_improp = True
        continue
    if line.startswith("MASS"):
        mass = True 
        continue
    if line.startswith("BOND"):
        bonds = True
        continue 
    if line.startswith("NONBOND") or line.startswith("MOD4"):
        i = 3
        while True:
            if cache[len(cache)-i] != "\n":
                dat = getVals(cache[len(cache)-i], "/")
                alias_dict[dat[0]] = dat
            else:
                break
            i += 1
        nonbond = True
        continue 
    if line == "\n" or line.strip(" ") == "\n":
        dihedrals_prop = False
        dihedrals_improp = False
        angles = False
        mass = False
        nonbond = False
        bonds = False
        continue
    if angles: 
        dat = parse_inputline_angle(line)
        angle_input_dict_parm[dat[0]+dat[1]+dat[2]] = dat[3:]
    if dihedrals_prop:
        dat = parse_inputline_dihedral_prop(line)
        names = dat[0:4]
        name = ""
        for a in names:
            if a != "X":
                name+=a
        if name not in dihedral_prop_input_dict_parm.keys():
            dihedral_prop_input_dict_parm[name] = [dat[4:]]
        else:
            dihedral_prop_input_dict_parm[name].append(dat[4:])
    
    if dihedrals_improp:
        dat = parse_inputline_dihedral_improp(line)

        names = dat[0:4]
        name = ""
        for a in names:
            if a != "X":
                name+=a
        if name not in dihedral_improp_input_dict_parm.keys():
            dihedral_improp_input_dict_parm[name] = [dat[4:]]
        else:
            dihedral_improp_input_dict_parm[name].append(dat[4:])
    if mass:
        dat = getVals(line, "/")[0:2]
        mass_nonbond_input_dict_parm[dat[0]] = [dat[1]]
    if nonbond: 
        dat = getVals(line, "/")[:3]
        if dat[0] in alias_dict.keys():
            for elem in alias_dict[dat[0]]:
                if elem in mass_nonbond_input_dict_parm.keys():
                    mass_nonbond_input_dict_parm[elem].append(convertR(dat[1]))
                    mass_nonbond_input_dict_parm[elem].append(convertEps(dat[2]))
                else:
                    print(f"Atom {elem} from nonbond parameters not found in mass directive")
        elif dat[0] in mass_nonbond_input_dict_parm.keys():
            mass_nonbond_input_dict_parm[dat[0]].append(convertR(dat[1]))
            mass_nonbond_input_dict_parm[dat[0]].append(convertEps(dat[2]))
        else:
            print(f"Atom {dat[0]} from nonbond parameters not found in mass directive")
    if bonds:
        dat = parse_inputline(line)
        bond_input_dict_parm[dat[0]+dat[1]] = dat[2:]
print("Secondary FF input file processed")

# Check atomtypes for matches

all_atomtypes_found = True
all_atomtypes_match = True
max_diff_sig = 0


not_found_atoms = []

for atomt in atomtype_list:
    not_found = False
    if atomt[0] in mass_nonbond_input_dict.keys() and len(mass_nonbond_input_dict[atomt[0]]) == 3:
        res = mass_nonbond_input_dict[atomt[0]]
    elif atomt[0] in mass_nonbond_input_dict_parm.keys() and len(mass_nonbond_input_dict_parm[atomt[0]]) == 3:
        res = mass_nonbond_input_dict_parm[atomt[0]]
    else:
        not_found_atoms.append(atomt[0])
        not_found = True
        all_atomtypes_found = False
        print(f"Atom {atomt[0]} not found in input!")
    if not not_found:
        max_diff_sig = max(max_diff_sig, abs(float(atomt[2])-float(res[1])))
        if float(atomt[1]) != float(res[0]) or float(atomt[3]) != float(res[2]):
            print(f"Atom {atomt[0]} has non-matching mass or epsillon!")
            all_atomtypes_match = False


clean_data = []
surviving_indices = []
current_index = 0
for moleculetype in data_raw:
    moleculetype_okay = True
    for atom in not_found_atoms:
        if atom in moleculetype["atoms"].values():
            moleculetype_okay = False
            print(f"Atom {atom} appears in a moleculetype, this moleculetype has been excluded!")
    if moleculetype_okay:
        clean_data.append(moleculetype)
        surviving_indices.append(current_index)
    current_index +=1

bond_list = []
angle_list = []
dihedral_list = []
for moleculetype in clean_data:
    for key in moleculetype.keys():
        if key == "bonds":
            bond_list.extend(moleculetype[key])
        if key == "angles":
            angle_list.extend(moleculetype[key])
        if key == "dihedrals":
            dihedral_list.extend(moleculetype[key])

# Check dihedrals for matches 

all_dihedrals_found = True
all_dihedrals_match = True
res = [-1,-1,-1,-1]
max_diff_dih = 0
for dihedral in dihedral_list:
    type_of_dih = dihedral[4]
    is_proper = (int(type_of_dih) != 4)
    perm = get_permutations_core(dihedral[0:4])
    perm_imp = get_permutations_4_wild(dihedral[0:4])
    found_prop = False
    for permu in perm:
        if not is_proper:
            break   
        if permu in dihedral_prop_input_dict.keys():
            cur = dihedral_prop_input_dict[permu]
            res = cur[findMatch(dihedral,cur)]
            found_prop = True 
            break
    found_improp = False
    if not found_prop:
        for permu in perm_imp:
            if is_proper:
                break
            if permu in dihedral_improp_input_dict.keys():
                cur = dihedral_improp_input_dict[permu]
                res = cur[findMatch(dihedral,cur)]
                found_improp = True
                break
    found_prop_2 = False 
    found_improp_2 = False
    if not found_improp and not found_prop:
        for permu in perm:
            if not is_proper:
                break
            if permu in dihedral_prop_input_dict_parm.keys():
                cur = dihedral_prop_input_dict_parm[permu]
                res = cur[findMatch(dihedral,cur)]
                found_prop_2 = True 
                break
        if not found_prop_2:
            for permu in perm_imp:
                if is_proper:
                    break
                if permu in dihedral_improp_input_dict_parm.keys():
                    cur = dihedral_improp_input_dict_parm[permu]
                    res = cur[findMatch(dihedral,cur)]
                    found_improp_2 = True
                    break
                    
    if not found_prop and not found_improp and not found_prop_2 and not found_improp_2:
        print(f"{dihedral} Dihedral not found anywhere!")
        all_dihedrals_found = False
    else:
        for i in range(0,4):
            if i == 1:
                max_diff_dih = max(max_diff_dih, abs(float(dihedral[5])-float(res[1])))
                dihedral[5] = res[1]
                continue
            if not float(dihedral[4+i]) == float(res[i]):
                print(f"Deviation {dihedral} and {cur} with {permu}")
                all_dihedrals_match = False
                break

# Check angles for matches

all_angles_found = True
all_angles_match = True 
max_diff = 0

for angle in angle_list: 
    perm = get_permutations_3(angle[0:3])
    not_in = True
    for permu in perm:
        if permu in angle_input_dict.keys():
            not_in = False
            if float(angle[5]) != float(angle_input_dict[permu][1]):
                print(f"Angle is off {angle}")
                all_angles_match = False
            max_diff = max(max_diff, abs(float(angle[4])-float(angle_input_dict[permu][0])))
            angle[4] = angle_input_dict[permu][0]
            break
    if not_in:
        for permu in perm: 
            if permu in angle_input_dict_parm.keys():
                not_in = False
                if float(angle[5]) != float(angle_input_dict_parm[permu][1]):
                    all_angles_match = False
                    print("f Angle is off {angle}")
                max_diff = max(max_diff, abs(float(angle[4])-float(angle_input_dict_parm[permu][0])))
                angle[4] = angle_input_dict_parm[permu][0]
                break
    if not_in:
        print(f"Angle {angle} not found!!")
        all_angles_found = False


# Check bonds for matches

all_bonds_found = True 
all_bonds_match = True

for bond in bond_list:
    current_a = bond[0]+bond[1]
    current_b = bond[1]+bond[0]
    if current_a in bond_input_dict.keys():
        if float(bond[3]) != float(bond_input_dict[current_a][0]) or float(bond[4]) != float(bond_input_dict[current_a][1]):
            all_bonds_match = False
            print(f"Bond {bond} is off!")
    elif current_b in bond_input_dict.keys():
       if float(bond[3]) != float(bond_input_dict[current_b][0]) or float(bond[4]) != float(bond_input_dict[current_b][1]):
           all_bonds_match = False
           print(f"Bond {bond} is off!")
    elif current_a in bond_input_dict_parm.keys():  
        if float(bond[3]) != float(bond_input_dict_parm[current_a][0]) or float(bond[4]) != float(bond_input_dict_parm[current_a][1]):
            all_bonds_match = False            
            print(f"Bond {bond} is off!")
    elif current_b in bond_input_dict_parm.keys():
        if float(bond[3]) != float(bond_input_dict_parm[current_b][0]) or float(bond[4]) != float(bond_input_dict_parm[current_b][1]):
            all_bonds_match = False
            print(f"Bond {bond} is off!")
    else:
        all_bonds_found = False
print("\n")
print("Results for finding all values:")
print("\n")

error = False 

if not all_atomtypes_found: 
    print("Some atomtypes were missing. Read the lines above. If any of these atomtypes do NOT belong to solvent or ions something is inconsistent with the provides forcefield, otherwise this is fine.")
else: 
    print("All atomtypes were found in input.")

if not all_bonds_found:
    error = True
    print("ERROR: Some bonds were missing from input.")
else:
    print("All bonds were found in input.")

if not all_angles_found:
    error = True
    print("ERROR: An angle was not found in the input.")
else:
    print("All angles were found in input.")

if not all_dihedrals_found:
    error = True
    print("ERROR: Some dihedrals were missing.")
else:
    print("All dihedrals were found in input.")

print("\n")
print("Results for comparing all values:")
print("\n")

if not all_atomtypes_match:
    error = True
    print("ERROR: Some atomtype values (except sigma) are off.")
else:
    print("All atomtype values (except sigma) match.")

if not all_bonds_match:
    error = True
    print("ERROR: Some bond values are off.")
else:
    print("All bond values match.")

if not all_angles_match:
    error = True
    print("ERROR: Some angle values (except angles) are off.")
else:
    print("All angle values (except angles) match.")

if not all_dihedrals_match:
    error = True
    print("ERROR: Some dihedral values (except angles) are off.")
else:
    print("All dihedral values (except angles) match.") 

print(f"Maximum deviation in angles of angles is {max_diff}")
print(f"Maximum deviation in angles of dihedrals is {max_diff_dih}.")
print(f"Maximum deviation in sigma values is {max_diff_sig}.")
print("\n")

if error:
    print("An error occured (see above), you should find out what causes is before it makes sense to fix the angles")
    exit(0)

fixed_top_name = get_file_no_suffix(top_file_name) + "_fixed.top"

fixed_top_file = open(fixed_top_name, "w")
top_file = open(top_file_name, "r")

angles = False
dihedrals = False

angle_count = 0
dihedral_count = 0


if max_diff > 0 or max_diff_dih > 0:
    molecule_index = -1
    for line in top_file:
        if not angles and not dihedrals:
            fixed_top_file.write(line)
        if line.startswith("[ moleculetype ]"):
            molecule_index += 1
        if line.startswith("[ angles ]") and molecule_index in surviving_indices:
            angles = True
            continue 
        if line.startswith("[ dihedrals ]") and molecule_index in surviving_indices:
            dihedrals = True
            continue 
        if (angles or dihedrals) and line.startswith(";"):
            fixed_top_file.write(line)
            continue
        if (angles or dihedrals) and (line == "\n" or line.strip(" ") == "\n"):
            fixed_top_file.write(line)
            angles = False
            dihedrals = False
            continue 
        if angles:
            replacement = str(float(angle_list[angle_count][4]))
            line = replace_angle(line, getVals(line,";")[4], replacement)
            fixed_top_file.write(line)
            angle_count += 1
        if dihedrals:
            replacement = str(float(dihedral_list[dihedral_count][5])) 
            line = replace_angle(line, getVals(line,";")[5], replacement)
            fixed_top_file.write(line)
            dihedral_count += 1

    print(f"Angles succesfully corrected in file {fixed_top_name}, sigma value deviations can not be corrected exactly, but are minor.")
else:
    print("All angles are exact. No correction needed!")

print("\n")
gro_fix_wished = input("Do you want to match coordinates of two .gro files (one generated from AMBER, one generated from GROMACS)? This is not generally not necessary (the slight differences are changed during MD anyway) but can be useful to compare energies of initial states. The two files must have the same atomnumbers and should not contain solvent or ions. [y/n]\n")
if gro_fix_wished == "y":
    gromacs_file_name = input("Provide the .gro file you have generated via GROMACS (with suffix) \n") 
    amber_file_name = input("Provide the .gro file you have generated via AMBER (with suffix) \n")
    out_file_name = get_file_no_suffix(amber_file_name) + "_fixed.gro"
    amber_file = open(amber_file_name, "r+")
    gromacs_file = open(gromacs_file_name, "r+")
    out_file = open(out_file_name, "w")
    start = True
    while True: 
        amber_line = amber_file.readline()
        gromacs_line = gromacs_file.readline()
        if start:
            out_file.write(amber_line)
            start = False
            continue
        if gromacs_line == "":
            out_file.write(amber_line)
        else:
            amber_values = {}
            gromacs_values = {}
            amber_split = amber_line.split(" ")
            gromacs_split = gromacs_line.split(" ")
            for i in range(0,len(amber_split)):
                if amber_split[i] == "":
                    continue 
                else: 
                    amber_values[i] = amber_split[i]
            for i in range(0, len(gromacs_split)):
                if gromacs_split[i] == "":
                    continue 
                else:
                    gromacs_values[i] = gromacs_split[i]
            if len(gromacs_values.keys()) == 3:
                out_file.write(amber_line)
                continue 
            cnt = 0
            for key_1,key_2 in zip(gromacs_values.keys(), amber_values.keys()):
                if cnt < 3:
                    cnt+=1 
                    continue 
                else:
                    amber_split[key_2] = gromacs_values[key_1]
            out_file.write(" ".join(amber_split))
        if amber_line == "":
            break



