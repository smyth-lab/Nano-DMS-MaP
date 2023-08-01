import xml.etree.ElementTree as ET
import numpy as np

###################### Generate csv and bp2seq files ###############


def read_in_xml(xml_file, sample, with_stdev = False):

    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    transcript_id = root[0].attrib["id"]
    length = root[0].attrib["length"]
    sequence = root[0][0].text.replace("\t", "").replace("\n", "")
    reactivity = np.array(root[0][1].text.replace("\t", "").replace("\n", "").split(",")).astype(float)
    if with_stdev:
        stdev = np.array(root[0][2].text.replace("\t", "").replace("\n", "").split(",")).astype(float)
        return {"sample" : sample, 
            "transcript_id" : transcript_id,
            "length" : length,
            "sequence" : sequence, 
            "reactivity" : reactivity,
            "stdev": stdev
           }
    else:
        return {"sample" : sample, 
            "transcript_id" : transcript_id,
            "length" : length,
            "sequence" : sequence, 
            "reactivity" : reactivity
           }
    
    
def convert_xml_to_bpseq(xml_file,outfile):

    tmp_data = read_in_xml(xml_file, "")

    reactivities = tmp_data["reactivity"]
    sequence = list(tmp_data["sequence"].replace("T", "U"))
    
    reactivities = np.nan_to_num(reactivities, nan=-1.0)
    with open(outfile, "w+") as out:
        for i in np.arange(1,1+reactivities.shape[0]):
            position = int(i)
            line = f"{position} {sequence[position-1]} e1 {reactivities[position-1]}\n"
            out.write(line)
            
            

def convert_xml_to_alt_formats(data_folder, sample, isoform, option_list = ["q22_eq10_ndni"], reactive_nt_list = ["AC", "ACT", "ACGT"], norm_option_list = [""]):    

    for option in option_list: #["q22_eq10_ndni", "q22_eq10", "default"]

        #specify reactive_nt again as in rfnorm
        for reactive_nt in reactive_nt_list: #["ACGT", "AC", "ACT", "G"]

            #include "_raw" if rfnorm was also run with Zubradt (4)
            for norm_option in norm_option_list: #["", "_raw"]

                xml_combined = f"{data_folder}/rfnorm/{sample}/{isoform}/{option}_{reactive_nt}{norm_option}/{isoform}.xml"
                if os.path.isfile(xml_combined):
                    file = read_in_xml(xml_combined, sample)
                    reactivities = file["reactivity"]
                    reactivities = np.nan_to_num(reactivities, nan = -1)
                    reactivity_file = xml_combined.replace(".xml", ".csv")
                    np.savetxt(reactivity_file, reactivities, fmt="%6f")
                    bp_file = xml_combined.replace(".xml", ".bp2seq")
                    convert_xml_to_bpseq(xml_combined, bp_file)

                    
                    

################## MaP File Generation #########################
#this optional step requires "raw", i.e. unnormalized reactivities not subtracted by control. Required rf-norm settings are Zubradt and "-r", i.e. raw output. In addition, specify 6 decimals (-d 6) during rf-combine
#

import numpy as np
def stderr(mutrate, coverage):
    return np.sqrt(mutrate) / np.sqrt(coverage)
    
def stderr_position(stderr_sample, stderr_untreated):
    return np.sqrt((stderr_sample)**2 + (stderr_untreated)**2)


def convert_xml_to_map(data_folder, sample, control, isoform, option_list = ["q22_eq10_ndni"], reactive_nt_list = ["AC", "ACT", "ACGT"]):
    
    coverage_file = f"{os.getcwd()}/data/bam/{sample}/{sample}_{isoform}_coverage.csv"
    coverage_data = pd.read_csv(coverage_file, sep="\t", names=["RNA", "position", "coverage"])

    control_coverage_file = f"{os.getcwd()}/data/bam/{control}/{control}_{isoform}_coverage.csv"
    control_coverage_data = pd.read_csv(control_coverage_file, sep="\t", names=["RNA", "position", "coverage"])
    positions = coverage_data["position"].values
    coverage = coverage_data["coverage"].values
    control_coverage = control_coverage_data["coverage"].values
    
    for option in option_list:

        for nts in reactive_nt_list:
            os.makedirs(f"{data_folder}/map_files/{nts}/{isoform}", exist_ok=True)
            xml_file = f"{os.getcwd()}/data/rfnorm/{sample}/{isoform}/{option}_{nts}_raw/{isoform}.xml"
            control_xml_file = f"{os.getcwd()}/data/rfnorm/{control}/{isoform}/{option}_{nts}_raw/{RNA}.xml"
            try:
                xml_data = read_in_xml(xml_file, sample)
            except:
                print("Could not find xml file", xml_file)
                continue
            try:
                control_xml_data = read_in_xml(control_xml_file, sample)
            except:
                print("Could not find xml file", control_xml_file)
                continue
            
            reactivity = xml_data["reactivity"]
            sequence = np.array(list(xml_data["sequence"]))
            
            control_reactivity = control_xml_data["reactivity"]
            
            stderr_sample = np.sqrt(reactivity) / np.sqrt(coverage)
            stderr_control = np.sqrt(control_reactivity) / np.sqrt(control_coverage)
            stderrs = np.sqrt(stderr_sample**2 + stderr_control**2)
            
            diff_reactivity = reactivity - control_reactivity
            
            tmp_df = pd.DataFrame(np.array([positions, np.nan_to_num(diff_reactivity, nan=-999), np.nan_to_num(stderrs, nan=0), sequence]).T, columns = ["position", "diff_raw_reactivity", "stderr", "nt"])
            tmp_df["diff_raw_reactivity"] = tmp_df["diff_raw_reactivity"].astype(float)
            tmp_df["stderr"] = tmp_df["stderr"].astype(float)
            outfile = f"{os.getcwd()}/data/rfnorm/{sample}/{RNA}/{option}_{nts}_raw/{RNA}.map"
            tmp_df.round(6).to_csv(outfile, sep="\t", float_format='%.6f', header=False, index=False)
