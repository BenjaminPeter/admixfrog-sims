import yaml
from collections import Iterable
from pprint import pprint
from itertools import product

#series = "S_RLEC"

#C = yaml.load(open("sims.yaml"))

def get_iter(C, kw, series):
    R = C[kw]['__default__'].copy()
    for k, v in R.items():
        if not isinstance(v, Iterable):
            R[k] = [v]
    for k, v in C['series'][kw][series].items():
        if not isinstance(v, Iterable):
            C['series'][kw][series][k] = [v]

    #pprint(R)
    #pprint(C['series'][kw][series])
    
    R.update(C['series'][kw][series])
    return R #product(*R.values())

def get_series_runs(R, series):
    for i, vals in enumerate(product(*R.values())):
        D = dict((name, v) for name, v in zip(R.keys(), vals))
        s_name = f'{series}_{i}'
        yield s_name, D


def update_all(C):
    for kw in C['series']:
        for series in C['series'][kw]:
            R = get_iter(C, kw, series)
            C[kw].update(dict(get_series_runs(R, series)))

    

#S = get_iter(C, 'sim', series)
#D = get_iter(C, 'demography', series)
#A = get_iter(C, 'admixfrog', series)
#R = get_iter(C, 'rle', series)


