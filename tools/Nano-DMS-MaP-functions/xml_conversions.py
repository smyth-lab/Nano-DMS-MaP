import xml.etree.ElementTree as ET
import numpy as np

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
            
            
#function to help predicting only part of an RNA isoform
def convert_xml_to_bpseq_trimmed(xml_file,outfile, length):

    tmp_data = read_in_xml(xml_file, "")

    reactivities = tmp_data["reactivity"][:length]
    sequence = list(tmp_data["sequence"].replace("T", "U"))[:length]
    
    reactivities = np.nan_to_num(reactivities, nan=-1.0)
    with open(outfile, "w+") as out:
        for i in np.arange(1,1+reactivities.shape[0]):
            position = int(i)
            line = f"{position} {sequence[position-1]} e1 {reactivities[position-1]}\n"
            out.write(line)
            