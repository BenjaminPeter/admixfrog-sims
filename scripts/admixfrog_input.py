import numpy as np
import pandas as pd
from scipy.stats import binom, poisson
def admixfrog_input(snps, rec, coverage, contamination,chrom='1', prefix='admixfrog',
                    ascertainment=None, all_deni=False):
    if 'pos' not in snps:
        snps['pos'] = snps.index.astype(int)
    snps.pos = snps.pos.astype(int)
    snps2 = snps.drop_duplicates(subset=['pos'])
    afr_cols = [col for col in snps2.columns if col.startswith("afr")] 
    asn_cols = [col for col in snps2.columns if col.startswith("asn")] 
    eur_cols = [col for col in snps2.columns if col.startswith("eur")] 
    cont_cols = [col for col in snps2.columns if col.startswith("cont")] 

    den_cols = [col for col in snps2.columns if col.startswith("den")]
    if all_deni:
        den_cols = den_cols[1:]
    else:
        den_cols = den_cols[3:]


    samples = asn_cols + eur_cols + den_cols
    n_samples = len(samples) // 2
    samples = np.array(samples).reshape(n_samples, 2).tolist()

    n_afr, n_asn, n_eur = len(afr_cols), len(asn_cols), len(eur_cols)
    n_den, n_cont = len(den_cols), len(cont_cols)

    D = dict()
    D['chrom'] = chrom
    D['pos'] = snps2.pos.astype(int)
    if rec is None:
        D['map'] = snps2.pos / 1e6
    else:
        D['map'] = np.interp(snps2.pos, rec.pos, rec.map)
    D['ref'] = 'A'
    D['alt'] = 'G'


    D["AFR_alt"] = np.sum(snps2[afr_cols], 1)
    D["ALT_alt"] = snps2.nea0 + snps2.nea1
    D["CHA_alt"] = snps2.nea2 + snps2.nea3
    D["VIN_alt"] = snps2.nea4 + snps2.nea5
    D["D11_alt"] = binom.rvs(1, p= (snps2.nea6 + snps2.den0) / 2.)
    D["DEN_alt"] = snps2.den1 + snps2.den2
    D["CONT_alt"] = np.sum(snps2[cont_cols], 1)

    D["ALT_ref"] = 2 - D['ALT_alt']
    D["CHA_ref"] = 2 - D['CHA_alt']
    D["VIN_ref"] = 2 - D['VIN_alt']
    D["DEN_ref"] = 2 - D['DEN_alt']
    D["D11_ref"] = 1 - D['D11_alt'] #pseudo-haploid
    D["AFR_ref"] = n_afr - D['AFR_alt']
    D["CONT_ref"] = n_cont - D['CONT_alt']

    D["PAN_alt"] = snps2.pan0 * 2
    D["PAN_ref"] = 2 - D['PAN_alt']

    D["NEA_ref"] = D['ALT_ref'] + D['VIN_ref']
    D["NEA_alt"] = D['ALT_alt'] + D['VIN_alt']
    D["D_ref"] = D['DEN_ref'] + D['D11_ref']
    D["D_alt"] = D['DEN_alt'] + D['D11_alt']


    if n_asn:
        D["ASN_alt"] = np.sum(snps2[asn_cols], 1)
        D["ASN_ref"] = n_asn - D['ASN_alt']
    ref = pd.DataFrame.from_dict(D)
    assert ref.shape[0] == snps2.shape[0]
    assert ref.drop_duplicates(subset=['chrom', 'pos']).shape[0] == ref.shape[0]

    if ascertainment == "archaicadmixture" :
        filter_ = (ref.AFR_alt + ref.ALT_ref == 0) | (ref.AFR_ref + ref.ALT_alt == 0)
        ref = ref[filter_]
        snps2 = snps2[filter_]
    elif ascertainment == "archaics":
        filter_ = (ref.NEA_alt + ref.DEN_alt + ref.AFR_alt != 0) & (ref.AFR_ref + ref.NEA_ref + ref.DEN_ref != 0)
        ref = ref[filter_]
        snps2 = snps2[filter_]
    elif ascertainment == "none":
        pass
    else:
        raise ValueError("ascertainment not known")

    #ref.to_csv(f"{prefix}.panel.xz", float_format="%.5f", index=False, compression="xz")
    ref.to_csv(f"{prefix}.panel.txt", float_format="%.5f", index=False)

    n_libs = len(coverage)
    assert len(contamination) == n_libs
    libs = [f"lib{i}" for i in range(n_libs)]

    for i, sample_ids in enumerate(samples):
        admixfrog_sample(sample_ids, 
                         ref=ref, 
                         snps=snps2,
                         coverage=coverage, 
                         contamination=contamination, 
                         libs=libs, 
                         name=f'{prefix}.{i}')


def admixfrog_ids(snps, prefix='admixfrog', all_deni=False):
    snps2 = snps
    afr_cols = [col for col in snps2.columns if col.startswith("afr")] 
    asn_cols = [col for col in snps2.columns if col.startswith("asn")] 
    eur_cols = [col for col in snps2.columns if col.startswith("eur")] 
    cont_cols = [col for col in snps2.columns if col.startswith("cont")] 

    den_cols = [col for col in snps2.columns if col.startswith("den")]
    if all_deni:
        den_cols = den_cols[1:]
    else:
        den_cols = den_cols[3:]


    samples = asn_cols + eur_cols + den_cols
    n_samples = len(samples) // 2
    samples = np.array(samples).reshape(n_samples, 2).tolist()

    df = pd.DataFrame([(i, s) for i, ss in enumerate(samples) for s in ss])
    df.columns = ['id', 'sample_id']                                       
    df.to_csv(prefix, index=False)

def admixfrog_sample(ids, ref, snps, coverage, contamination, libs, name):
    S = []
    for cov, cont, lib in zip(coverage, contamination, libs):
        print(f'Sample{name}\tLib:{lib}\tCov:{cov}\tcont:{cont}', end="\t")
        data = ref[['chrom', 'pos']].copy()
        data['true_alt'] = np.sum(snps[ids],1)
        data['true_ref'] = 2 - data['true_alt']
        data['lib'] = lib

        
        cov_real = poisson.rvs(cov * (1. - cont), size=data.shape[0])
        cov_cont = poisson.rvs(cov * cont, size=data.shape[0])
        p = data['true_alt'] / (data['true_ref'] + data['true_alt'])
        p_cont = ref.CONT_alt / (ref.CONT_ref + ref.CONT_alt)
        data['ralt'] = binom.rvs(cov_real, p)
        data['rref'] = cov_real - data['ralt']

        data['calt'] = binom.rvs(cov_cont, p_cont)
        data['cref'] = cov_cont - data['calt']

        data['talt'] = data.ralt + data.calt
        data['tref'] = data.rref + data.cref

        print(f"alt:\t{lib}\t{np.mean(data['ralt']):.3f}\t{np.mean(data['calt']):.3f}\t{np.mean(data['talt']):.3f}")
        print(f"ref:\t{lib}\t{np.mean(data['rref']):.3f}\t{np.mean(data['cref']):.3f}\t{np.mean(data['tref']):.3f}")

        data = data[data.tref+data.talt>0]
        S.append(data)
    data = pd.concat(S).sort_values(['chrom', 'pos', 'lib'])
    data.to_csv(f"{name}.sample.txt", float_format="%.5f",
                index=False
                )

debug_str = ""\
"--neand 55000,25,0,0.03 "\
"--eur-ages 45000 45000 40000 40000 30000 30000 "\
"--seq-len 5_000_000 "\
"--stats true_neand "\
"--introgression neand-eur "\
"--nafrA 40 "\
"--nafrB 20 "\
"--nasn 40 "\
"--snps "\
"--admixfrog "\
"--coverage 1 1 1 1 1 1 1 1 1 1 "\
"--contamination 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 "\
"--output-prefix tmp/cont".split()
#"--vcf --vcf-ploidy 1 "\
