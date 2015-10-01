deq = 0.335
S = float(660)
B = float(230)
num_helices = 28


# n is the number of basepairs installed in region with
#   a default of nref base pairs
# e.g. nref = 105
n = [0 for i in range(num_helices)]

# delta is the distance in nm from an arbitrary axis
# Assume 128 bp average length of helix
num_helices = 28
delta = [0 for i in range(num_helices)]
for helix_num in range(4):
	delta[helix_num] = 5
	n[helix_num] = 128 - 16
for helix_num in range(4, 10):
	delta[helix_num] = 3
	n[helix_num] = 128 - 10
for helix_num in range(10, 14):
	delta[helix_num] = 1
	n[helix_num] = 128 - 3
for helix_num in range(14, 18):
	delta[helix_num] = -1
	n[helix_num] = 128 + 3
for helix_num in range(18, 24):
	delta[helix_num] = -3
	n[helix_num] = 128 + 10
for helix_num in range(24, 28):
	delta[helix_num] = -5
	n[helix_num] = 128 + 16
for i in range(num_helices):
	delta[i] *= -1.125


total_n = 0
for i in range(num_helices):
	total_n += n[i]
average_n = float(total_n)/num_helices
print "average n is", average_n, "\n"


# Calculate delta_div_n_sum, one_div_n_sum, and beta values
delta_div_n_sum = 0
one_div_n_sum = 0
for i in range(num_helices):
	delta_div_n_sum += delta[i]/n[i]
	one_div_n_sum += float(1)/n[i]

beta = [0 for i in range(num_helices)]
for i in range(num_helices):
	beta[i] = (delta[i]-delta_div_n_sum/one_div_n_sum)/n[i]

# Calculate theta
beta_n_sum = 0
beta_delta_sum = 0
for i in range(num_helices):
	beta_n_sum += beta[i]*n[i]
	beta_delta_sum += beta[i]*delta[i]
numerator = deq*beta_n_sum
denominator_term_0 = beta_delta_sum
denominator_term_1 = (B/S)*(one_div_n_sum)
denominator = denominator_term_0 + denominator_term_1
theta = numerator/denominator

theta_in_degrees = theta*180/3.1416
print "theta is", theta_in_degrees,"degrees or ", theta, "radians"

# Calculate d
d = [0 for i in range(num_helices)]
for i in range(num_helices):
	factor_0 = theta*(delta[i] - delta_div_n_sum/one_div_n_sum)
	factor_1 = deq*num_helices/one_div_n_sum
	d[i] = (factor_0 + factor_1)/n[i]

for i in range(num_helices):
	print "helix", i, "\tlength per bp is", d[i]

print
print "delta_div_n_sum is", delta_div_n_sum
print "one_div_n_sum is", one_div_n_sum
print "beta_n_sum is", beta_n_sum
print "beta_delta_sum is", beta_delta_sum
print "numerator is", numerator
print "denominator_term_0 is", denominator_term_0
print "denominator_term_1 is", denominator_term_1
print "denominator is", denominator





