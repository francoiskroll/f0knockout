# Sample sizes simulations

Refers to Figure 7B.

___

Please cite the preprint if you use some of the data or code! <br />
https://www.biorxiv.org/content/10.1101/2020.06.04.133462v3

___

Files <br />
* *trpa1b_ifsomewts_SampleSize_100fish10sim.RData*
* *csnk1db_ifsomewts_SampleSize_100fish10sim.RData*

can be found at the Zenodo version of this repository (too heavy for GitHub):

https://doi.org/10.5281/zenodo.3900611

## Simulations for *trpa1b*: mustard oil response

Files in folder **trpa1b_simulation**.

Script `trpa1b_ifsomewts_SampleSize.R` uses *dpxresults.csv* (same as in **fig5a_trpa1b_mustard** folder) and runs the simulations.
`simulation_grid.pdf` is an example of the density plots for one simulation (i.e. 101 plots: from 0/100 KO in F0 group until 100/100 KO in F0 group).

As the simulations take a while to run (at least on my laptop), `trpa1b_ifsomewts_SampleSize.R` saves the Environment (*trpa1b_ifsomewts_SampleSize_100fish10sim.RData*).
Script `trpa1b_ifsomewts_SampleSize_fromRData.R` takes the RData as input and produces the barplot in Figure 7B left (*f0_trpa1b_samplesize.pdf*).

## Simulations for *csnk1db*: circadian period

Files in folder **csnk1db_simulation**.

Script `csnk1db_ifsomewts_SampleSize.R` uses *220719_periodBoth.csv* (same as in **fig5b_csnk1db_per3luc** folder) and runs the simulations.
`simulation_grid.pdf` is an example of the density plots for one simulation (i.e. 101 plots: from 0/100 KO in F0 group until 100/100 KO in F0 group).

As the simulations take a while to run (at least on my laptop), `csnk1db_ifsomewts_SampleSize.R` saves the Environment (*csnk1db_ifsomewts_SampleSize_100fish10sim.RData*).
Script `csnk1db_ifsomewts_SampleSize_fromRData.R` takes the RData as input and produces the barplot in Figure 7B right (*f0_csnk1db_samplesize.pdf*).

---

Feel free to get in touch for questions

  * [![alt text][1.2]][1] [@francois_kroll](https://twitter.com/francois_kroll)

  * :email: francois@kroll.be

<!-- icons with padding -->
[1.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)

<!-- icons without padding -->
[1.2]: http://i.imgur.com/wWzX9uB.png (twitter icon without padding)

<!-- links to your social media accounts -->
[1]: https://twitter.com/francois_kroll
