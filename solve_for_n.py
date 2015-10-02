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

num_potential_strands_range = [4, 20]
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

min_theta_error_candidate = {"theta_error": 100}

def print_candidate(num_segments, n_diff, num_potential_strands, n_diff_period):
    print('\n')

    theta = calc_bend_angle(num_segments*7, n_diff)
    theta_error = (360-theta)/3.6

    print str(num_potential_strands) + ' strands per well,',
    print str(num_segments) + ' segments,',
    print "{0:.2f}".format(theta_error) + '% theta error,',
    print str(num_segments*7) + ' bp long,',
    print str(n_diff) + " insertions/deletions ,",
    print str(n_diff/num_segments*n_diff_period) + ' insertions/deletions per',
    if n_diff_period == 1:
        print "segment"
    else:
        print str(n_diff_period) + ' segments'

    #calc_bend_angle(num_segments*7, n_diff, True, False, False)

    if abs(theta_error) < min_theta_error_candidate["theta_error"]:
        min_theta_error_candidate["theta_error"] = abs(theta_error)
        min_theta_error_candidate["num_segments"] = num_segments
        min_theta_error_candidate["n_diff"] = n_diff
        min_theta_error_candidate["n_diff_period"] = n_diff_period
        min_theta_error_candidate["num_potential_strands"] = num_potential_strands


tolerance = 5

def factors(num):
    facts =  set(reduce(list.__add__, ([i, num//i] for i in range(1, int(num**0.5) + 1) if num % i == 0)))
    return list(facts)

for num_potential_strands in range(num_potential_strands_range[0], num_potential_strands_range[1]+1):
    min_num_segments = int(math.ceil(num_segments_range[0]/num_potential_strands)*num_potential_strands)
    for num_segments in range(min_num_segments, num_segments_range[1]+1, num_potential_strands):

        n = num_segments*num_bp_per_segment
        n_diff = solve_for_n_diff(n)

        num_segments_factors = factors(num_segments)
        for i, val in enumerate(num_segments_factors):
            remainder = n_diff/val - round(n_diff/val)
            if abs(remainder*val) < tolerance and (round(n_diff/val) in num_segments_factors):
                factor1 = val
                factor2 = round(n_diff/val)
                n_diff = factor1*factor2
                print_candidate(num_segments, n_diff, num_potential_strands, num_segments/max(factor1, factor2))


print "\n\n"
print "candidate with min theta error:"
print_candidate(min_theta_error_candidate["num_segments"], min_theta_error_candidate["n_diff"],
                min_theta_error_candidate["num_potential_strands"], min_theta_error_candidate["n_diff_period"])
print calc_bend_angle(min_theta_error_candidate["num_segments"]*num_bp_per_segment, min_theta_error_candidate["n_diff"], True, False, True)
print "\n"
