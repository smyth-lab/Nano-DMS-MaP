from sklearn.metrics import roc_auc_score
import subprocess

def generate_db_from_eterna(eterna_outfile, sample, outdir, reactivity_file, mode="sample", start = 0, end = "full", mask_PBS = False, mask_DIS = False, primer_5 = 0, primer_3 = 0):
    db_files = []
    
    tmp_data = read_in_xml(reactivity_file, "")
    reactivity = tmp_data["reactivity"]
    sequence = tmp_data["sequence"]
    if end == "full":
        end = len(reactivity)
    else:
        if primer_3 > 0:
            print("You both wanted to trim sequence from the 3' end and mask reactivities due to a primer binding site. Are you sure?")

    if mask_PBS:
        reactivity[184:199] = 1.0
    
    if mask_DIS:
        reactivity[256:261] = 1.0
    
    if (start > 0) & (primer_5 > 0):
        print("You both wanted to trim sequence from the 5' end and mask reactivities due to a primer binding site. Are you sure?")
    
    reactivity[:primer_5] = -1.0
    if primer_3 > 0:
        reactivity[-primer_3:] = -1.0
    reactivity = reactivity[start:end]
    sequence = sequence[start:end]

    mask = np.isnan(reactivity)
    reactivity = reactivity[~mask]

    if mode == "sample":
        with open(eterna_outfile, "r") as infile:
            for i, line in enumerate(infile):
                structure = np.array(list(line.strip()))

                rocauc = np.round(roc_auc_score(np.array(list(structure))[~mask] == ".", reactivity),3)

                with open(f"{outdir}/{sample}_{i}.db", "w") as outfile:
                    outfile.write(f">{sample}_{i}_ROCAUC{rocauc}\n")
                    outfile.write(f"{sequence}\n")
                    outfile.write(f"{''.join(structure)}")
                db_files.append(f"{outdir}/{sample}_{i}.db")
        return db_files
    elif mode == "predict":
        with open(eterna_outfile, "r") as infile:
            _, _, _, _, structure = infile.readlines()
            structure = np.array(list(structure.strip()))
            rocauc = np.round(roc_auc_score(np.array(list(structure))[~mask] == ".", reactivity),3)

            with open(eterna_outfile.replace(".eterna", ".db"), "w") as outfile:
                outfile.write(f">{sample}_{tmp_data['transcript_id']}_predict_ROCAUC{rocauc}\n")
                outfile.write(f"{sequence}\n")
                outfile.write(f"{''.join(structure)}")
            db_files.append(eterna_outfile.replace(".eterna", ".db"))
        return db_files
    
def generate_varna(db_file, reactivity_file, sample, varna_outprefix, start = 0, end = "full", primer_5 = 0, primer_3 = 0):
    reactivities = read_in_xml(reactivity_file, "")["reactivity"]
    reactivities = np.nan_to_num(reactivities, nan=-1.0)
    
    if end == "full":
        end = len(reactivities)
    
    reactivities = reactivities[start:end]
    reactivities[:primer_5] = -1.0
    if primer_3 > 0:
        reactivities[-primer_3:] = -1.0
    #print(reactivities)
    
    with open(db_file, "r") as infile:
        title = infile.readline().strip()
        sequence = infile.readline().strip().replace("T", "U")
        structure = infile.readline().strip()
    
    for algorithm in ["radiate", "line"]:
        varna_outfile = f"{varna_outprefix}_{algorithm}.varna"
        
        colormap = '-1.0:#888888;0.0:#0000FF;0.5:#FFFFFF;1.0:#FF0000'
        command = f'java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -algorithm {algorithm} -sequenceDBN "{sequence}" -structureDBN "{structure}" -o {varna_outfile} -colorMap "{";".join(reactivities.astype(str))}" -colorMapStyle "{colormap}" -title "{title}" -flat True'
        #print(command)
        subprocess.run(command)
    return