import numpy as np


def smart_gen_array(func, a, b, sing=[], n=100, eps=1e-17):
    """
    Evaluate function func between a and b and replace singularities
    with np.nan. Every subinterval between the singularities is filled
    with n points and we evaluate up to eps around the singularity.
    """
    # make intervals between singularities
    intervals = []
    sorted_sing = sorted(sing)

    if len(sorted_sing) == 0:
        intervals.append([a, b])
    else:
        for i, s in enumerate(sorted_sing):
            if s < a:
                continue

            # first interval
            if s > a and len(intervals)==0:
                intervals.append([a, s])

            # last interval
            elif s > b:
                intervals.append([sorted_sing[i-1], b])

            # intermediate intervals
            else:
                intervals.append([sorted_sing[i-1], s])

    if (intervals[-1][1] < b):
        intervals.append([intervals[-1][1], b])

    # solve function within the intervals and stitch with np.nan
    # on the singularities

    xx = []
    yy = []

    for interval in intervals:
        a, b = interval

        # safety margin
        a += eps
        b -= eps

        x = np.linspace(a, b, n)
        y = func(x)

        x = np.append(x, b)
        y = np.append(y, np.nan)

        xx = np.append(xx, x)
        yy = np.append(yy, y)

    return np.array(xx), np.array(yy)


def laser_intensity(laser_power, beam_waist):
    """
    Laser intensity in Watts/cm^2 for laser power in W and and beam waist in m
    """
    return 2/np.pi * laser_power/(beam_waist/1e-2)**2
