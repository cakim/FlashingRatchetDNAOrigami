import sys

#for now we're assuming attach only to 3' end, but that can change

# Read in json file
input_file = file('caDNAno_sequences_raw.json', 'r')
sequence_data = eval(input_file.read())
input_file.close()
input_file = file('sticky_strands.json', 'r')
sticky_strands = eval(input_file.read())
input_file.close()

output = {}
output["scaffold"] = {
    "length": str(len(sequence_data["scaffold"])) + " of p7308",
    "seq": sequence_data["scaffold"]
}

sticky_staples = []
for staple_data in sequence_data["staples"]:
    seq = staple_data["seq"]
    seq += sticky_strands[staple_data["color"]]
    sticky_staples.append(seq)

output["staples"] = sticky_staples

import json
text_file = open("sequences_final.json", "w")
text_file.write(json.dumps(output, indent=4))
text_file.close()