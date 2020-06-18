# *scn1lab* stable knockout (behaviour tracking)

Refers to Figure 6B,C

___

Please cite the preprint if you use some of the data or code! <br />
https://www.biorxiv.org/content/10.1101/2020.06.04.133462v3

___

Files <br />
* *140218_34_DATA.mat*
* *190813_06_DATA.mat*
* *190813_07_DATA.mat*

can be found at the Zenodo version of this repository (too heavy for GitHub):

https://doi.org/10.5281/zenodo.3900611

Main script is MATLAB script `Vp_Analyse_fingerprint.m`, which takes three mat files as input (in that order):
* *140218_34_DATA.mat* (same as in folder **scn1lab_stable**)
* *190813_06_DATA.mat* (same as in folder **scn1lab_f0**)
* *190813_07_DATA.mat* (same as in folder **scn1lab_f0**)

It generates files:
* *exptags.csv*
* *grptags.csv*
* *fingerprint.csv*
* *eucdists.csv*

*exptags.csv*; *grptags.csv*; *fingerprint.csv* are input to script `fingerprint.R`, which generates the fingerprint plot and the correlation heatmap in Figure 6B (*f0_fingerprint.pdf*, *f0_corrheat.pdf*).

*exptags.csv*; *grptags.csv*; *eucdists.csv* are input to script `euclideanscatter.R`, which generates the Euclidean distance scatterplot in Figure 6C (*f0_eucdist.pdf*). It also performs statistical tests.

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
