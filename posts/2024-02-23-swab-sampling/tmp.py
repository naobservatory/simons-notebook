#!/usr/bin/env python3

FRAC_SICK=0.01
SIMULATIONS=10000

import math
import random
from collections import defaultdict

# Reads and Covid Reads columns from Lu et al 2021 Table S1
# pbpaste | tr ',' '_' | while read a b ; do echo "($a, $b)," ; done
LU_DATA = [
    # reads, covid reads
    (42_084_470, 33_275_395),
    (27_707_635, 9_770_670),
    (45_408_518, 3_351_784),
    (56_421_340, 464_316),
    (46_154_092, 617_330),
    (68_776_916, 2_318_220),
    (38_439_824, 831_304),
    (16_691_298, 101_478),
    (31_118_799, 32_108),
    (43_935_756, 35_525),
    (38_421_258, 11_456),
    (85_901_880, 10_997),
    (24_343_922, 17_482),
    (43_891_267, 8_328),
    (31_491_439, 150),
]

def probabilistic_round(x):
    return int(math.floor(x + random.random()))

def simulate_once(n_swabs):
    n_sick = probabilistic_round(FRAC_SICK*n_swabs)
    if not n_sick:
        return 0
    ra = 0
    for _ in range(n_sick):
        reads, covid_reads = random.choice(LU_DATA)
        ra += covid_reads / reads
    return ra / n_swabs

def simulate_many(n_simulations, n_swabs):
    RAs = []
    for _ in range(n_simulations):
        RAs.append(simulate_once(n_swabs))
    RAs.sort()

    return [RAs[probabilistic_round(len(RAs)*percentile/100)]
            for percentile in range(100)]

data = defaultdict(list)
for n_swabs in [50, 100, 200, 400, 800, 1000000]:
    for percentile, ra in enumerate(simulate_many(SIMULATIONS, n_swabs)):
        data[percentile].append(ra)
    
for percentile, ras in sorted(data.items()):
    print (percentile, *ras)
