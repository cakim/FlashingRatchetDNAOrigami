# coding=utf-8


def calc_bend_angle(nRef, nDiff, degrees, verbose):

    deq = 0.335 # unperturbed length-per-base-pair in B-form DNA

    # one that includes twist-stretch coupling (S=1100 pN, B=230 pN· nm, D=-90pN· nm)
    # and another set (S=660 pN, B=230 pN· nm, D=0 pN· nm) excluding twist-stretch coupling and with S and B
    # calculated via a Young’s modulus derived from a DNA persistence length of 50 nm and our experimentally
    # observed effective helix diameter
    S = float(660)
    B = float(230)


    num_helices = 6
    helix_diam = 2.25


    # n is the number of basepairs installed in region with
    #   a default of nref base pairs
    # e.g. nref = 105
    n = [0 for i in range(num_helices)]

    # delta is the distance in nm from an arbitrary axis
    # n is number of bp in each helix
    delta = [0 for i in range(num_helices)]
    for helix_num in range(2):
        delta[helix_num] = -helix_diam
        n[helix_num] = nDiff - nDiff
    for helix_num in range(2, 4):
        delta[helix_num] = 0
        n[helix_num] = nDiff
    for helix_num in range(4, 6):
        delta[helix_num] = helix_diam
        n[helix_num] = nDiff + nDiff

    print(n)
    total_n = 0
    for i in range(num_helices):
        total_n += n[i]
    average_n = float(total_n)/num_helices
    if (verbose):
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
    if (verbose):
        print "theta is", theta_in_degrees,"degrees or ", theta, "radians"

    # Calculate d
    d = [0 for i in range(num_helices)]
    for i in range(num_helices):
        factor_0 = theta*(delta[i] - delta_div_n_sum/one_div_n_sum)
        factor_1 = deq*num_helices/one_div_n_sum
        d[i] = (factor_0 + factor_1)/n[i]

    if (verbose):
        for i in range(num_helices):
            print "helix", i, "\tlength per bp is", d[i]

        print "delta_div_n_sum is", delta_div_n_sum
        print "one_div_n_sum is", one_div_n_sum
        print "beta_n_sum is", beta_n_sum
        print "beta_delta_sum is", beta_delta_sum
        print "numerator is", numerator
        print "denominator_term_0 is", denominator_term_0
        print "denominator_term_1 is", denominator_term_1
        print "denominator is", denominator

    if (degrees):
        return theta_in_degrees
    return theta



