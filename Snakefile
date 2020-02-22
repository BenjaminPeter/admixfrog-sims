from pprint import pprint
import pandas as pd
from scripts.expand_series import update_all
from scripts.admixfrog_input import admixfrog_input, admixfrog_ids
configfile : "config/sims.yaml"
include: 'plots.snake'

#localrules: run_sim, sample_table, merge_panel, merge_sample, merge_true, classify_frags, hapmap_rec,

C=config

N_SAMPLES =  8
N_REPS = 20 
CHROMS = range(1,21)
#CHROMS = range(1, 23)
#CHROMS = range(15, 23)

update_all(config)
#pprint(config['admixfrog'])

"""
basic outline
1. generate simulation

2.1 (Truth)
    a. get all true fragments

3.1 (Data) 
    a. generate admixfrog input files
    b. generate admixfrog panels
    c. merge panels / input chromosome

4. run admixfrog
    a. run admixfrog
    b. call fragments

5. (evaluation)
    a. compare fragments with truth
    b. compare contamination

"""


def get_rec_true(wc):
    sim = config['sim']['__default__'].copy()
    sim.update(config['sim'][wc.sim])
    map_ = sim['rec']
    return f"recs/hapmap/{map_}/chr{wc.chrom}.rec.gz"
rule run_sim:
    input:
        rec=get_rec_true
    priority : 1
    group : "sim"
    benchmark:
        "benchmarks/run_sim/{sim}/{chrom}/{demo}.{rep}_snps.tsv",
    output:
        data="sims/{sim}/{chrom}/{demo}.{rep}_snps.tsv.gz",
        frags="sims/{sim}/{chrom}/{demo}.{rep}_all_haplotypes.csv.gz",
    run:
        """1. get demography and simulation settings from config"""
        demo = config['demography']['__default__'].copy()
        demo.update(config['demography'][wildcards.demo])
        sim = config['sim']['__default__'].copy()
        sim.update(config['sim'][wildcards.sim])

        """2. set up gene flow events"""
        gf_nea_eur = ' --gf-nea-eur {t_gf_neaeur},{d_gf_neaeur},{m_gf_neaeur},{m_gf_neaeur_rev} '.format(**demo)
        gf_nea_den = ' --gf-nea-den {t_gf_neaden},{d_gf_neaden},{m_gf_neaden},{m_gf_neaden_rev} '.format(**demo)
        gf_den_eur = ' --gf-den-eur {t_gf_deneur},{d_gf_deneur},{m_gf_deneur},{m_gf_deneur_rev} '.format(**demo)

        """3. set up samples """
        """   3.1reference panel sizes for african, asians and oceanians"""
        n = ' --n-afr {n_afr} --n-asn {n_asn} --n-oce {n_oce} '.format(**demo)

        """   3.2 european (early modern human test samples )  """
        if 'a_eur' in demo:
            n += ' --a-eur ' + " ".join(f"{i} {i}" for i in demo['a_eur'])
        else:
            n += ' --n-eur 0 ' 

        """   3.3 Denisovan (archaic test samples )  """
        if 'a_den' in demo:
            n += ' --a-den ' + " ".join(f"{i} {i}" for i in demo['a_den'])

        exe = '../demography/demography_deni3.py'
        pars = '--writesnps  --introgression nea-eur eur-nea nea-den den-nea den-eur eur-den' 
        opt = ' --output-prefix sims/{wildcards.sim}/{wildcards.chrom}/{wildcards.demo}.{wildcards.rep}'
        chrom = " --chrom {wildcards.chrom} "
        #rec = " --rec {input.rec} "
        rec = ""

        seq = '--seq-len {seq_len}'.format(**sim)

        s = " ".join((exe ,gf_nea_eur, gf_nea_den, gf_den_eur, rec, n, chrom, pars, opt, seq))
        shell(s)

rule sample_table:
    input:
        snps="sims/{sim}/%s/{demo}.0_snps.tsv.gz" % CHROMS[0],
    output:
        sample_table="sims/{sim}/{demo}.idtbl",
    run:
        snps = pd.read_csv(input.snps, sep="\t")
        sim = config['sim']['__default__'].copy()
        sim.update(config['sim'][wildcards.sim])

        admixfrog_ids(snps, 
            prefix=output.sample_table)
        

def get_rec_obs(wc):
    sim = config['coverage']['__default__'].copy()
    sim.update(config['coverage'][wc.cov])
    map_ = sim['rec_obs']
    return f"recs/hapmap/{map_}/chr{wc.chrom}.rec.gz"
rule admixfrog_input:
    """create admixfrog input and panel file for 1 run 1 chromosome"""
    input:
        snps="sims/{sim}/{chrom}/{demo}.{rep}_snps.tsv.gz",
        rec=get_rec_obs,
    resources:
        io=1
    priority : 0
    group : "sim"
    output:
        panel=temp("sims/{sim}/{chrom}/{demo}.{cov}.{rep}.panel.txt"),
        samples=temp(expand("sims/{{sim}}/{{chrom}}/{{demo}}.{{cov}}.{{rep}}.{id}.sample.txt", 
            id=range(N_SAMPLES)))
    run:
        cov = config['coverage']['__default__'].copy()
        cov.update(config['coverage'][wildcards.cov])

        snps = pd.read_csv(input.snps, sep="\t")
        rec = pd.read_csv(input.rec, sep="\t")

        coverage = [c[0] for c in cov['coverage']]
        contamination = [c[1] for c in cov['coverage']]
        asc = cov['ascertainment']

        admixfrog_input(snps=snps, rec=rec, chrom=wildcards.chrom, 
            coverage=coverage,
            contamination=contamination,
            ascertainment = asc,
            prefix=".".join(output.panel.split(".")[:-2]))

rule merge_panel:
    input:
        expand("sims/{{sim}}/{chrom}/{{demo}}.{{cov}}.{{rep}}.panel.txt", chrom=CHROMS)
    priority : 100
    group : "sim"
    output:
        "infiles/{sim}/{demo}.{cov}.{rep}.panel.xz"
    shell:
        "head -qn1 {input[0]} | xz -c > {output} && "
        "tail -qn+2 {input} | xz -c >> {output} "

rule merge_sample:
    input:
        expand("sims/{{sim}}/{chrom}/{{demo}}.{{cov}}.{{rep}}.{{id}}.sample.txt", chrom=CHROMS)
    priority : 1020
    group : "sim"
    output:
        "infiles/{sim}/{demo}.{cov}.{rep}.{id}.sample.xz"
    shell:
        "head -qn1 {input[0]} | xz -c > {output} && "
        "tail -qn+2 {input} | xz -c >> {output} "

def get_recs(wc):
    pars = config['sim']['__default__'].copy()
    pars.update(config['sim'][wc.sim])
    return expand("recs/hapmap/{map}/chr{chrom}.rec.gz",
                 map=pars['rec'],
                 chrom = CHROMS)
rule merge_true:
    input:
        f=expand("sims/{{sim}}/{chrom}/{{demo}}.{{rep}}_all_haplotypes.csv.gz", chrom=CHROMS),
        _script = 'scripts/merge_true.R'
        #rec_files = get_recs,
    priority : 10
    group : "sim"
    output:
        frags="sims/{sim}/{demo}.{rep}_all_haplotypes.xz",
    script: "scripts/merge_true.R"

rule merge_true_reps:
    input:
        f=expand("sims/{{sim}}/{{demo}}.{rep}_all_haplotypes.xz",
            rep=range(N_REPS)),
        idtbl ='sims/{sim}/{demo}.idtbl'
    output:
        frags="sims/true/{sim}/{demo}_all_haplotypes_merged.xz"
    script: 'scripts/merge_true2.R'

rule run_admixfrog:
    input:
        sample="infiles/{sims}/{demo}.{cov}.{rep}.{id}.sample.xz",
        ref="infiles/{sims}/{demo}.{cov}.{rep}.panel.xz",
    params:
        freq_f = 3,
        freq_c = 1,
        ancestral = "PAN",
        max_iter = 200,
        n_post_rep = 100,
        ll_tol = 0.01,
        error = 0.01,
    output:
        f=expand("admixfrog/{{pars}}/{{sims}}/{{demo}}.{{cov}}.{{rep}}.{{id}}.{ext}.xz",
            ext=['bin', 'cont', 'res', 'res2'])
    group: 'frags'
    run:
        wc = wildcards
        pars = config['admixfrog']['__default__'].copy()
        pars.update(config['admixfrog'][wc.pars])
        state_str = " ".join(pars['state_ids'])

        s = "admixfrog --infile {input.sample} --ref {input.ref} -o  "
        s += ' admixfrog/{wc.pars}/{wc.sims}/{wc.demo}.{wc.cov}.{wc.rep}.{wc.id}'
        s += f" --states {state_str} "
        s += f" --cont-id {pars['cont_id']} "
        s += " --ll-tol {params.ll_tol} "
        s += f" --bin-size {pars['bin_size']} "
        s += " --est-F --est-tau --freq-F {params.freq_f} "
        s += " --freq-contamination {params.freq_c} "
        s += " --e0 {params.error} "
        s += " --est-error "
        s += " --ancestral {params.ancestral} "
        s += " --max-iter {params.max_iter}"
        s += " --n-post-replicates {params.n_post_rep}"
        s += " --no-snp --no-rle > /dev/null"
        
        print(s)
        shell(s)

rule admixfrog_frags:
    input:
        bins="admixfrog/{pars}/{sims}/{demo}.{cov}.{rep}.{id}.bin.xz",
    output:
        rle="rle/{pars}/{rle}/{sims}/{demo}.{cov}.{rep}.{id}.rle.xz",
    group: 'frags'
    run:
        pars = config['rle']['__default__'].copy()
        pars.update(config['rle'][wildcards.rle])
        penalty = pars['run_penalty']

        s = "admixfrog-rle --in {input.bins} --out {output.rle} "
        s += " --run-penalty {penalty} "
        shell(s)
        

rule classify_frags:
    input:
        estfile="rle/{pars}/{rle}/{sim}/{demo}.{cov}.{rep}.{id}.rle.xz",
        truefile="sims/{sim}/{demo}.{rep}_all_haplotypes.xz",
        sample_table="sims/{sim}/{demo}.idtbl",    
    output:
        frags = "admixfrog/{pars}/{sim}/{rle}/{demo}.{cov}.{rep}.{id}.frags"
    script: "scripts/compare_frags.R"

def _all_simfrags(wc):
    demo = [x for x in C['demography'] if x.startswith(f"{wc.demo_series}_")]
    sim = [x for x in C['sim'] if x.startswith(f"{wc.sim_series}_")]
    pars = [x for x in C['admixfrog'] if x.startswith(f"{wc.af_series}_")]
    cov = [x for x in C['coverage'] if x.startswith(f"{wc.cov_series}_")]
    frags=expand("admixfrog/{pars}/{sim}/{demo}.{cov}.{rep}.{id}.res.xz",
        pars=pars,
        sim=sim,
        demo=demo,
        cov=cov,
        rep=range(N_REPS),
        id=range(N_SAMPLES))
    if len(frags) == 0:
        raise ValueError("no input files")
    return frags
rule merge_simfrags_series:
    input:
        frags=_all_simfrags,
        script_="scripts/merge_simruns.R",
    output:
        frags="series/simfrags/{af_series}/{sim_series}/{demo_series}/{cov_series}.xz"
    script: "scripts/merge_simruns.R"

def _all_series(wc):
    rle = [x for x in C['rle'] if x.startswith(f"{wc.rle_series}_")]
    demo = [x for x in C['demography'] if x.startswith(f"{wc.demo_series}_")]
    sim = [x for x in C['sim'] if x.startswith(f"{wc.sim_series}_")]
    pars = [x for x in C['admixfrog'] if x.startswith(f"{wc.af_series}_")]
    cov = [x for x in C['coverage'] if x.startswith(f"{wc.cov_series}_")]
    frags = expand("admixfrog/{pars}/{sim}/{rle}/{demo}.{cov}.{rep}.{id}.frags",
        pars=pars,
        sim=sim,
        demo=demo,
        rle=rle,
        cov=cov,
        rep=range(N_REPS),
        id=range(N_SAMPLES))
    if len(frags) == 0:
        raise ValueError("no input files")
    return frags
rule classify_series:
    input:
        frags = _all_series,
        _script = "scripts/merge_frags.R"
    output:
        frags="series/{af_series}/{sim_series}/{rle_series}/{demo_series}/{cov_series}.class.xz"
    script: "scripts/merge_frags.R"

def _all_est_series(wc):
    rle = [x for x in C['rle'] if x.startswith(f"{wc.rle_series}_")]
    demo = [x for x in C['demography'] if x.startswith(f"{wc.demo_series}_")]
    sim = [x for x in C['sim'] if x.startswith(f"{wc.sim_series}_")]
    pars = [x for x in C['admixfrog'] if x.startswith(f"{wc.af_series}_")]
    cov = [x for x in C['coverage'] if x.startswith(f"{wc.cov_series}_")]
    frags = expand("rle/{pars}/{rle}/{sim}/{demo}.{cov}.{rep}.{id}.rle.xz",
        pars=pars,
        sim=sim,
        demo=demo,
        rle=rle,
        cov=cov,
        rep=range(N_REPS),
        id=range(N_SAMPLES))
    if len(frags) == 0:
        raise ValueError("no input files")
    return frags
rule merge_est_series:
    """all estimated fragments for a series"""
    input:
        frags = _all_est_series,
        _script = "scripts/merge_est_frags.R"
    output:
        frags="series/{af_series}/{sim_series}/{rle_series}/{demo_series}/{cov_series}.est.xz"
    script: "scripts/merge_est_frags.R"

def _all_cont(wc):
    demo = [x for x in C['demography'] if x.startswith(f"{wc.demo_series}_")]
    sim = [x for x in C['sim'] if x.startswith(f"{wc.sim_series}_")]
    pars = [x for x in C['admixfrog'] if x.startswith(f"{wc.af_series}_")]
    cov = [x for x in C['coverage'] if x.startswith(f"{wc.cov_series}_")]
    frags = expand("admixfrog/{pars}/{sim}/{demo}.{cov}.{rep}.{id}.cont.xz",
        pars=pars,
        sim=sim,
        demo=demo,
        cov=cov,
        rep=range(N_REPS),
        id=range(N_SAMPLES))
    if len(frags) == 0:
        raise ValueError("no input files")
    return frags
rule classify_cont:
    input:
        frags = _all_cont,
        _script = "scripts/merge_cont.R"
    output:
        frags="series/{af_series}/{sim_series}/{demo_series}/{cov_series}.rds"
    script: "scripts/merge_cont.R"


rule hapmap_rec:
    input:
        rec="recs/maps_b37/maps_chr.{chrom}",
        script_= "scripts/haprec.R"
    output:
        rec="recs/hapmap/{map}/chr{chrom}.rec.gz"
    script: "scripts/haprec.R"



def series_cfg(wc):
    A = config['set'][wc.name]['admixfrog']
    S = config['set'][wc.name]['sim']
    R = config['set'][wc.name]['rle']
    D = config['set'][wc.name]['demography']
    C = config['set'][wc.name]['coverage']
    
    s1 = f'series/{A}/{S}/{R}/{D}/{C}.class.xz'
    s2 = f'series/{A}/{S}/{R}/{D}/{C}.est.xz'
    s3 = f'series/simfrags/{A}/{S}/{D}/{C}.xz'
    s4 = f'series/{A}/{S}/{D}/{C}.rds'
    return s1, s2, s3, s4
rule run_series_cfg:
    input:
        series_cfg
    output:
        touch("seriesx/{name}")
