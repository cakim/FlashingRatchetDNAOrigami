from solve_for_num_insertions_deletions import solve_for_n_diff
from calculate_bend_angle import calc_bend_angle
import math

#constraints to satisfy:
#n_diff (number of insertions/deletions across the total track length) is integer
#n is a multiple of 7bp (num_segments is integer)
#n_diff/num_segments <= 1, no more than one deletion insertion per segment
#num_segments = n/7 is a multiple of the # of different strands in the potential well (between 4 and 7 to keep things reasonable)
#n is less than 700bp - we don't want to design a huge, huge thing
#num_segments/n_diff is close to integer - distribute deletions/insertions evenly

num_potential_strands_range = [4, 9]
num_segments_max = 100#n = 700bp

num_bp_per_segment = 7

#calculate minimum_num_segments (n_Diff/num_segments = 1, 1 deletion/insertion in every segment of tensioned/compressed helices)
def calculate_num_segments_range():

    for i in range(15, num_segments_max+1):
        n_diff = solve_for_n_diff(num_bp_per_segment*i)
        if n_diff/i <= 1:
            return [i, num_segments_max]

    print "problem computing min num segments, abort"
    return None

num_segments_range = calculate_num_segments_range()


min_theta_error_candidate = {"theta_erorr": 100}

def print_candidate(num_segments, n_diff, num_potential_strands):
    print('\n')

    num_seg_per_del_insert = round(num_segments/n_diff)
    final_n_diff = num_segments/num_seg_per_del_insert
    theta = calc_bend_angle(num_segments*7, final_n_diff)
    theta_error = (360-theta)/3.6

    print(str(num_potential_strands) + ' strands per well, ' + str(num_segments) + ' segments, ' +
          "{0:.2f}".format(theta_error) + '% theta error, ' + str(num_segments*7) + ' bplong, ' +
          str(num_segments/final_n_diff) + ' insertions/deletions per segment')

    calc_bend_angle(num_segments*7, n_diff, True, False, True)

    if abs(theta_error) < min_theta_error_candidate["theta_erorr"]:
        min_theta_error_candidate["theta_erorr"] = abs(theta_error)
        min_theta_error_candidate["num_segments"] = num_segments
        min_theta_error_candidate["n_diff"] = final_n_diff
        min_theta_error_candidate["num_potential_strands"] = num_potential_strands


tolerance = 0.1

for num_potential_strands in range(num_potential_strands_range[0], num_potential_strands_range[1]+1):
    min_num_segments = int(math.ceil(num_segments_range[0]/num_potential_strands)*num_potential_strands)
    for num_segments in range(min_num_segments, num_segments_range[1]+1, num_potential_strands):
        n = num_segments*num_bp_per_segment
        n_diff = solve_for_n_diff(n)
        remainder = num_segments/n_diff - round(num_segments/n_diff)
        if abs(remainder) < tolerance:
            print_candidate(num_segments, n_diff, num_potential_strands)

print "\n\n"
print "candidate with min theta error:"
print_candidate(min_theta_error_candidate["num_segments"], min_theta_error_candidate["n_diff"],
                min_theta_error_candidate["num_potential_strands"])
print "\n"

#n_diff = solve_for_n_diff(n, 360)
#
#
#theta = calc_bend_angle(n, n_diff, True, True)
#print n_diff
#print n_diff*7/n