# QuickLook
quicklook light curve generator for TESS targets

![img](tests/TOI1150_s55_pdcsap_sc.png)

## Use case
* Given target name, run periodogram on TESS lightcurve (if it exists) to estimate the stellar rotation period and the orbital period of a potential companion i.e. planet brown, dwarf, or star

## Installation
* Install editable version within a conda environment
```bash

$ conda activate my_env
(my_env) $ pip install -r requirements.txt
(my_env) $ pip install -e .
```

## Run
```bash
(my_env) $ ql
usage: ql [-h] [-name NAME] [-sec SECTOR] [-lc {pdcsap,sap}]
          [-p {spoc,tess-spoc,cdips,pathos,qlp,tglc}] [-e EXPTIME]
          [-fm FLATTEN_METHOD] [-pm {gls,ls,bls}] [-wl WINDOW_LENGTH]
          [-ec EDGE_CUTOFF] [--sigma_clip_raw SIGMA_CLIP_RAW SIGMA_CLIP_RAW]
          [--sigma_clip_flat SIGMA_CLIP_FLAT SIGMA_CLIP_FLAT]
          [-plims PERIOD_LIMITS PERIOD_LIMITS] [-s] [-o OUTDIR] [-v]
          [--overwrite] [-img]
          [--survey {dss1,poss2ukstu_red,poss2ukstu_ir,poss2ukstu_blue,poss1_blue,poss1_red,all,quickv,phase2_gsc2,phase2_gsc1}]
          [-em EPHEM_MASK EPHEM_MASK EPHEM_MASK]

create quick look image of TESS data

options:
  -h, --help            show this help message and exit
  -name NAME            target name
  -sec SECTOR, --sector SECTOR
                        TESS sector
  -lc {pdcsap,sap}, --fluxtype {pdcsap,sap}
                        type of lightcurve
  -p {spoc,tess-spoc,cdips,pathos,qlp,tglc}, --pipeline {spoc,tess-spoc,cdips,pathos,qlp,tglc}
                        lightcurve produced from which pipeline
  -e EXPTIME, --exptime EXPTIME
                        exposure time (default is whatever is found in last
                        sector)
  -fm FLATTEN_METHOD, --flatten_method FLATTEN_METHOD
                        wotan flatten method (default=biweight)
  -pm {gls,ls,bls}, --pg_method {gls,ls,bls}
                        periodogran method (default=gls)
  -wl WINDOW_LENGTH, --window_length WINDOW_LENGTH
                        flatten method window length (default=0.5 days)
  -ec EDGE_CUTOFF, --edge_cutoff EDGE_CUTOFF
                        cut each edges (default=0.1 days)
  --sigma_clip_raw SIGMA_CLIP_RAW SIGMA_CLIP_RAW
                        (sigma_lo,sigma_hi) for outlier rejection after
                        flattening lc
  --sigma_clip_flat SIGMA_CLIP_FLAT SIGMA_CLIP_FLAT
                        (sigma_lo,sigma_hi) for outlier rejection after
                        flattening lc
  -plims PERIOD_LIMITS PERIOD_LIMITS, --period_limits PERIOD_LIMITS PERIOD_LIMITS
                        period limits in TLS search; default=(0.5, baseline/2)
                        d
  -s, --save            save figure and tls
  -o OUTDIR, --outdir OUTDIR
                        output directory
  -v, --verbose         show details
  --overwrite           overwrite files
  -img, --use_archival_image
                        plot gaia sources on archival image instead of tpf
  --survey {dss1,poss2ukstu_red,poss2ukstu_ir,poss2ukstu_blue,poss1_blue,poss1_red,all,quickv,phase2_gsc2,phase2_gsc1}
                        archival image survey name if using img option
  -em EPHEM_MASK EPHEM_MASK EPHEM_MASK, --ephem_mask EPHEM_MASK EPHEM_MASK EPHEM_MASK
                        mask ephemeris given period and t0
```

## Tools
[pre-commit](https://pre-commit.com/)
