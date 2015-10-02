from calculate_bend_angle import calc_bend_angle


tolerance = 0.1
n = 500
target_theta = 120


nDiff = 0
nDiff_max = n
nDiff_min = 0
theta = 0


##binary search across nDiff to find right theta for length
#while abs(theta-target_theta) > tolerance:
#
#    if theta < target_theta:
#        nDiff_min = nDiff
#    else:
#        nDiff_max = nDiff
#
#    nDiff = (nDiff_max+nDiff_min)/2.0
theta = calc_bend_angle(250, 10, True, False)


print theta
print nDiff