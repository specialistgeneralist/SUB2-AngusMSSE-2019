A Statistical Timetable for the sub-2 hour Marathon
===================================================

*(To appear in MSSE, 2019)*

*Author: Simon D. Angus*

*Dept. of Economics, Monash University*

*Melbourne, Australia.*

## Requirements

You will need:
* **MATLAB** (tested with version r2018b); and
* **Statistics Toolbox** (provides *fitnlm*, *predict*, and *fitgev*).

## To replicate the paper

Open a **MATLAB** prompt, navigate to the working folder, and to run for example, the _male_ version of the code run:
```
main('wr_male.csv', [0.5 0.2 0.10 0.04 0.02 0.01], 1)
```
with arguments
* `'wr_male.csv'`: use the male CSV input file;
* `[0.5 0.2 0.10 0.04 0.02 0.01]`, the `ALPHA` values to use for the run;
* `1`: sets `ISMALE=1`.

To run the _female_ version, one runs equivalently:
```
main('wr_female.csv', [0.5 0.2 0.10 0.04 0.02 0.01], 0)
```

When run, *main* will produce two _gendergap_ files for later comparison with *gender_gap*. For instance, if `ISMALE=1` the two files will be: `out_expected_male.csv`, and `out_1in10_male.csv` . After running both male and female variants, one can then conduct a long-run gender gap analysis as follows:
```
OUTNAME = 'gg_1in10.csv';
gender_gap('out_1in10_male.csv','out_1in10_female.csv',OUTNAME)
```

## Help?

Please get in touch at `simon.angus@monash.edu`.
