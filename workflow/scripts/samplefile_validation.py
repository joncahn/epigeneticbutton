import re

def validations_patterns(data_type):
    if data_type == "RNAseq":
        return re.compile(r"^RNAseq$")
    elif data_type == "RAMPAGE":
        return re.compile(r"^RAMPAGE$")
    elif data_type == "sRNA":
        return re.compile(r"^(sRNA|smallRNA|shRNA)$")
    elif data_type == "mC":
        return re.compile(r"^(mC|WGBS|ONT|Pico|EMseq)$")
    elif data_type.startswith("TF_"):
        return re.compile(r"^(IP|IPb|Input)$")
    elif data_type.startswith("ChIP"):
        return re.compile(r"^(?!.*\s)(?!.*__)(?!.*').*$")

def validate_sample_type(row):
    pattern = validations_patterns(row.data_type)
    return bool(pattern and pattern.fullmatch(row.sample_type))

def validate_SRA(row):
    id = str(row.seq_id)
    path = str(row.fastq_path)
    if id.startswith("SRR"):
        return path == "SRA"
    if path == "SRA":
        return id.startswith("SRR")
    return True

def name(row):
    return f"{row.data_type}_{row.line}_{row.tissue}_{row.sample_type}"

def assign_chip_input(row, tab):
    dtype = row.data_type
    stype = row.sample_type
    if dtype.startswith("TF") or dtype.startswith("ChIP"): 
        if stype != "Input":
            match = tab[
                (tab["data_type"]==row.data_type) &
                (tab["line"]==row.line) &
                (tab["tissue"]==row.tissue) &
                (tab["sample_type"]=="Input") &
                (tab["ref_genome"]==row.ref_genome)
            ]
            if match.empty:
                return "Input"

        if stype == "Input":
            match = tab[
                (tab["data_type"]==row.data_type) &
                (tab["line"]==row.line) &
                (tab["tissue"]==row.tissue) &
                (tab["sample_type"] != "Input") &
                (tab["ref_genome"]==row.ref_genome)
            ]
            if match.empty:
                return "Sample"          
    return True

def check_table(tab):
    err = 0
    dup = tab[tab.duplicated(
        subset=["data_type","line","tissue","sample_type","replicate","ref_genome"], 
        keep=False
    )]
    
    if not dup.empty:
        for _, r in dup.iterrows():
            print(f'[X] Duplicated rows: {name(r)}', flush=True)
            err = 1

    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        if not validate_sample_type(row):
            print(f'[X] Row #{i} {name(row)}: sample_type does not match data_type', flush=True)
            err = 1
        if not validate_SRA(row):
            print(f'[X] Row #{i} {name(row)}: fastq_path should be "SRA" for SRR IDs or a directory otherwise', flush=True)
            err = 1
        result = assign_chip_input(row, tab)
        if result == "Input":
            print(f'[X] Row #{i} {name(row)}: missing a corresponding Input sample', flush=True)
            err = 1
        elif result == "Sample":
            print(f'[X] Row #{i} {name(row)}: no sample depends on this Input', flush=True)
            err = 1

    if err == 0:
        print("[OK] Samplefile is correct!", flush=True)
    else:
        raise ValueError("[X] Validation failed — please fix your samplefile and rerun.")

