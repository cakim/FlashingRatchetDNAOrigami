from calculate_bend_angle import calc_bend_angle

#constraints to satisfy:
#nDiff (number of insertions/deletions across the total track length) is integer


def solve_for_n_diff(n, target_theta=360):
    
    tolerance = 0.01

    n_diff = 0 # num insertions deletions
    last_n_diff = -1
    n_diff_max = n
    n_diff_min = 0
    theta = 0

    def tie_break(n_diff1, n_diff2):
        theta1 = calc_bend_angle(n, n_diff1)
        theta2 = calc_bend_angle(n, n_diff2)

        if abs(theta1-target_theta) < abs(theta2-target_theta):
            return n_diff1
        return n_diff2

    ##binary search across integer n_diff to find right theta for length
    while abs(theta-target_theta) > tolerance:
    
        if theta < target_theta:
            n_diff_min = n_diff
        else:
            n_diff_max = n_diff
    
        n_diff = round((n_diff_max+n_diff_min)/2.0)
        if last_n_diff == n_diff:
            if n_diff-1 < n_diff_min:
                if n_diff+1 > n_diff_max:
                    break
                else:
                    n_diff += 1
            else:
                n_diff -= 1
    
        if n_diff == n_diff_max or n_diff == n_diff_min:# we've already tried this one
            n_diff = tie_break(n_diff, last_n_diff)
            break
    
        theta = calc_bend_angle(n, n_diff)
        last_n_diff = n_diff

    return n_diff

