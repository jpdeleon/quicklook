# QuickLook
`quicklook` is a Python program that runs a simple pipeline to search for transit signal in TESS (and Kepler soon) light curves. This program can be run in a jupyter notebook (see [example](https://github.com/jpdeleon/quicklook/tree/main/notebook)) or from the terminal using the `ql` script.

## Use case
Given target name, run periodogram on a TESS or Kepler lightcurve (if it exists) to estimate the stellar rotation period and the orbital period of a potential companion i.e. planet, brown dwarf, or star.
Although `quicklook` is optimized to find transiting exoplanets, this tool can also find eclipsing binaries and many other periodic signals.

## Try it on Google colab

<a href="https://colab.research.google.com/github/jpdeleon/quicklook/blob/main/notebook/examples.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>


## Installation
Create a conda environment called, say `my_env`, and install there an editable version of `quicklook`
```bash
$ conda create -n my_env python=3.10
$ conda activate my_env
(my_env) $ python -m python -m pip install -r https://raw.githubusercontent.com/jpdeleon/quicklook/main/requirements.txt
(my_env) $ python -m pip install -e git+https://github.com/jpdeleon/quicklook.git#egg=quicklook
```

If you want to run `quicklook` in a notebook, you also need to install jupyter
```
(my_env) $ python -m pip install jupyterlab
```

## Command line script
```bash
(my_env) $ ql
usage: ql [-h] [-name NAME] [-sec SECTOR] [-lc {pdcsap,sap}] [-p {spoc,tess-spoc,tasoc,cdips,pathos,qlp,tglc}] [-e EXPTIME] [-fm FLATTEN_METHOD] [-pm {gls,ls,bls}] [-wl WINDOW_LENGTH] [-ec EDGE_CUTOFF]
          [--sigma_clip_raw SIGMA_CLIP_RAW SIGMA_CLIP_RAW] [--sigma_clip_flat SIGMA_CLIP_FLAT SIGMA_CLIP_FLAT] [-plims PERIOD_LIMITS PERIOD_LIMITS] [-s] [-o OUTDIR] [-v] [--overwrite] [-img]
          [--survey {dss1,poss2ukstu_red,poss2ukstu_ir,poss2ukstu_blue,poss1_blue,poss1_red,all,quickv,phase2_gsc2,phase2_gsc1}] [-em EPHEM_MASK EPHEM_MASK EPHEM_MASK]

create quick look image of TESS data

options:
  -h, --help            show this help message and exit
  -name NAME            target name
  -sec SECTOR, --sector SECTOR
                        TESS sector (default=-1 (last available sector))
  -lc {pdcsap,sap}, --fluxtype {pdcsap,sap}
                        type of lightcurve
  -p {spoc,tess-spoc,tasoc,cdips,pathos,qlp,tglc}, --pipeline {spoc,tess-spoc,tasoc,cdips,pathos,qlp,tglc}
                        lightcurve produced from which pipeline (default=SPOC)
  -e EXPTIME, --exptime EXPTIME
                        exposure time (default is whatever is found in last sector)
  -fm FLATTEN_METHOD, --flatten_method FLATTEN_METHOD
                        wotan flatten method (default=biweight)
  -pm {gls,ls,bls}, --pg_method {gls,ls,bls}
                        periodogran method (default=gls)
  -wl WINDOW_LENGTH, --window_length WINDOW_LENGTH
                        flatten method window length (default=0.5 days)
  -ec EDGE_CUTOFF, --edge_cutoff EDGE_CUTOFF
                        cut each edges (default=0.1 days)
  --sigma_clip_raw SIGMA_CLIP_RAW SIGMA_CLIP_RAW
                        (sigma_lo,sigma_hi) for outlier rejection of raw lc before flattening/detrending
  --sigma_clip_flat SIGMA_CLIP_FLAT SIGMA_CLIP_FLAT
                        (sigma_lo,sigma_hi) for outlier rejection of flattened/detrended lc
  -plims PERIOD_LIMITS PERIOD_LIMITS, --period_limits PERIOD_LIMITS PERIOD_LIMITS
                        period limits in TLS search; default=(0.5, baseline/2) d
  -s, --save            save figure and tls
  -o OUTDIR, --outdir OUTDIR
                        output directory
  -v, --verbose         show details
  --overwrite           overwrite files
  -img, --use_archival_image
                        plot gaia sources on archival image instead of tpf
  --survey {dss1,poss2ukstu_red,poss2ukstu_ir,poss2ukstu_blue,poss1_blue,poss1_red,all,quickv,phase2_gsc2,phase2_gsc1}
                        archival image survey name if using img option
  -em, --use_ephem_mask
                        mask transits using TFOP ephemeris if available (default=False)
```

## Examples

1. Run `quicklook` on TOI 1150.01's most recent TESS lightcurve

```shell
(my_env) $ ql -name TOI-1150 -img
```
![img](tests/TOI1150_s55_pdcsap_sc.png)

The figure above shows 9 panels. Let's break them down.
* top row
  - left (panel 1): raw lightcurve (blue marker) and trend (red line)
  - middle (panel 2): [Lomb-Scargle periodogram](https://docs.astropy.org/en/stable/timeseries/lombscargle.html) used to estimate the star's rotation period; this is useful to find active and variable stars
  - right (panel 3): raw lightcurve phase-folded at the computed peak of Lomb-Scargle periodogram (corresponding to the stellar rotation period) from panel 1;
* middle row
  - left (panel 4): flattened lightcurve and detected transits (determined from the TLS periodogram on panel 5)
  - middle (panel 5): periodogram using the [transit least squares](https://ui.adsabs.harvard.edu/abs/2019A%26A...623A..39H/abstract) (TLS) algorithm
  - right (panel 6): TESS aperture (blue polygon) and annotated Gaia sources (orange and red markers) overlaid on archival [DSS](https://archive.stsci.edu/cgi-bin/dss_form) image centered on the target; this is useful to see if there are potential sources of the signal other than the target
* bottom row
  - left (panel 7): phase-folded lightcurve at the derived peak of TLS periodogram (corresponding to the orbital period); odd (red markers) and even transits (blue markers) and best-fit transit model (black line) are also shown
  - middle (panel 8): phase-folded lightcurve zoomed at phase=0.5 to check for a secondary eclipse which is a strong indicator of a self-luminous companion such as an eclipsing binary or a high-albedo brown dwarf; the computed transit depth (dashed line) is shown for reference
  - right (panel 9): summary info about the star and potential companion (e.g. planet candidate)

Try changing the parameters:
```shell
(my_env) $ ql -name TIC52368076 -v -s
(my_env) $ ql -name TOI-125.01 -v -s -p qlp #specific pipeline
(my_env) $ ql -name TOI-125.01 -v -s -sec 2 #specific TESS sector
(my_env) $ ql -name TOI-125.01 -v -s -fm cosine #specific function to detrend baseline
(my_env) $ ql -name TOI-125.01 -v -s -plims 1 5 #limit search between 1-5 days
```

## Advanced usage

If you would like to run `ql` on a list of TIC IDs (saved as `tic_ids.txt`), then you can make a batch script named `run_ql_given_tic.batch`. The output files containing the plots (*.png) and periodogram results (*_tls.h5) will be saved in `tic_dir` directory:

```shell
(my_env) $ cat tic_ids.txt | while read tic; do echo ql -name $tic -s -o tic_dir; done > run_ql_given_tic.batch
```

To test the Nth line of the batch script,

```shell
(my_env) $ cat run_ql_given_tic.batch | sed -n Np | sh
```

To run all the lines in parallel using [GNU parallel](https://www.gnu.org/software/parallel/) with N cores,

```shell
(my_env) $ cat run_ql_given_tic.batch | parallel -j N
```

After the batch script is done running, we can rank `ql` output in terms of Signal Detection Efficiency (SDE, See [Hippke et al. 2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...623A..39H/abstract)) using `read_tls` script:

```shell
(my_env) $ read_tls tic_dir
```
