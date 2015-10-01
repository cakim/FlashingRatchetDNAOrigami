import sys

def comp_seq_FN(raw_sequence):
    uppercase = {'a':'A', 'A':'A', 'c':'C', 'C':'C', 'g':'G', 'G':'G', 't':'T', 'T':'T'}
    complement = {'a':'T', 'A':'T', 'c':'G', 'C':'G', 'g':'C', 'G':'C', 't':'A', 'T':'A'}
    antisense_seq = ''
    for letter in raw_sequence:
	if letter in uppercase:
            antisense_seq = complement[letter] + antisense_seq
    return antisense_seq


# This program will read in a json file and a sequence file, and assign the sequence
# to the scaffold, then generate staple strand sequences, along with annotation
# of colors, beginning base, end base, length


# Read in json file
input_file = file('test.json', 'r')
all_file = eval(input_file.read())
input_file.close()

# Read in m13 sequence (p7308)
input_file = file('p7308.txt', 'r')
scaf_seq = input_file.read()
input_file.close()

vstrands = all_file['vstrands']

num_vstrands = len(vstrands)
num_helices = num_vstrands
vsn = {}
for vstrand_num in range(num_vstrands):
	vsn[vstrands[vstrand_num]['num']] = vstrand_num
vsn[-1] = -1


# Assume that there is just one scaffold 5prime end
# Find that 5prime end
null_bp = [-1, -1]
num_bases = len(vstrands[0]['scaf'])
for helix_num in range(num_helices):
	for base_num in range(num_bases):
		bpp = vstrands[vsn[helix_num]]['scaf'][base_num]
		if (bpp[:2] == null_bp) and (bpp[2:] != null_bp):
			start_bp = [helix_num, base_num]


# Assign m13 sequence to the scaffold path
# Initialize data structure
scaf_base_seq_dc = {}
for helix_num in range(num_helices):
	scaf_base_seq_dc[helix_num] = ['.' for i in range(num_bases)]

# Fill in sequence
counter = 0
[helix_num, base_num] = start_bp
base_length = 1 + vstrands[vsn[helix_num]]['loop'][base_num]
base_length += vstrands[vsn[helix_num]]['skip'][base_num]
scaf_base_seq_dc[helix_num][base_num] = scaf_seq[counter:counter + base_length]
counter += base_length
while vstrands[vsn[helix_num]]['scaf'][base_num][2:] != null_bp:
	[helix_num, base_num] = vstrands[vsn[helix_num]]['scaf'][base_num][2:]
	base_length = 1 + vstrands[vsn[helix_num]]['loop'][base_num]
	base_length += vstrands[vsn[helix_num]]['skip'][base_num]
	scaf_base_seq_dc[helix_num][base_num] = scaf_seq[counter:counter + base_length]
	counter += base_length
print scaf_base_seq_dc

# Generate list of staple strand paths
# Search for all 5prime ends of stap paths
# For each of these, parse to the end of that path
stap_path_ra = []
for start_helix_num in range(num_helices):
	for start_base_num in range(num_bases):
		bpp = vstrands[vsn[start_helix_num]]['stap'][start_base_num]
		if (bpp[:2] == null_bp) and (bpp[2:] != null_bp):
			sub_ra = []
			[helix_num, base_num] = [start_helix_num, start_base_num]
			sub_ra.append([helix_num, base_num])
			while vstrands[vsn[helix_num]]['stap'][base_num][2:] != null_bp:
				[helix_num, base_num] = vstrands[vsn[helix_num]]['stap'][base_num][2:]
				sub_ra.append([helix_num, base_num])
			stap_path_ra.append(sub_ra)


stap_color_dc  = {}
stap_color_dc[13369344] = 'red'
stap_color_dc[16204552] = 'red orange'
stap_color_dc[16225054] = 'light orange'
stap_color_dc[11184640] = 'olive'
stap_color_dc[5749504]  = 'light green'
stap_color_dc[29184]    = 'dark green'
stap_color_dc[243362]   = 'cyan'
stap_color_dc[1507550]  = 'blue'
stap_color_dc[7536862]  = 'purple'
stap_color_dc[12060012] = 'magenta'
stap_color_dc[3355443]  = 'dark gray'
stap_color_dc[8947848]  = 'light gray'


# Initialize a color data structure
stap_color_2d_dc = {}
for helix_num in range(num_helices):
	stap_color_2d_dc[helix_num] = [-1 for i in range(num_bases)]
# Assign colors
for helix_num in range(num_helices):
	for [base_num, color_value] in vstrands[vsn[helix_num]]['stap_colors']:
		stap_color_2d_dc[helix_num][base_num] = stap_color_dc[color_value]
	

# Generate staple strand sequences, annotate with color, start base, end base
for sub_ra in stap_path_ra:
	[helix_num, base_num] = sub_ra[0]
	stap_color = stap_color_2d_dc[helix_num][base_num]
	start_base = sub_ra[0]
	end_base = sub_ra[-1]
	seq = ''
	for [helix_num, base_num] in sub_ra:
		seq += comp_seq_FN(scaf_base_seq_dc[helix_num][base_num])
	annotation_string  = str(stap_color) + ', ' + 'start ' + str(start_base) + ', '
	annotation_string += 'end ' + str(end_base) + ', length ' + str(len(seq))
	print seq + '\t\t' + annotation_string



# Color annotation


			










