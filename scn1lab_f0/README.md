# *scn1lab* (behaviour tracking)

Refers to Figure 6A and Figure 6—figure supplement 1.

___

Please cite the preprint if you use some of the data or code! <br />
https://doi.org/10.1101/2020.06.04.133462

___

See Materials and Methods/Behavioural data analysis.

*190813_06_07.zip*: the thousands of xls files after raw file Export from ZebraLab (see *trpa1b* mustard oil assay README).

*080819_06genotypeMap.xlsx*: layout of 96-well plate in box 6. <br />
*080819_07genotypeMap.xlsx*: layout of 96-well plate in box 7. <br />
Converted to *190819_06_genotype.txt* and *190819_07_genotype.txt* with `genotypeGenerator.R`.

*190813_06_DATA.txt* + the thousands of xls files + *190819_06_genotype.txt* are used as input to `VpExtract.m` (in folder **behaviour_matlabscripts**/**Vp_Extract**). It generates MATLAB data file *190813_06_DATA.mat*

*190813_07_DATA.txt* + the thousands of xls files + *190819_07_genotype.txt* are used as input to `VpExtract.m` (in folder **behaviour_matlabscripts**/**Vp_Extract**). It generates MATLAB data file *190813_07_DATA.mat*

Note: Script `VpExtract.m` was originally forked from https://github.com/ghoshm/Structure_Paper. The version included has some small modifications.

*190813_06_DATA.mat* and *190813_07_DATA.mat* are used together as input to *VpAnalyse.m* (in folder **behaviour_matlabscripts**/**Vp_Analyse**). It generates different plots not included in the paper and performs statistical tests (see *190813_both_anova.pdf*).

About statistics: in F0 paper, quote

> During the day, both F0 knockouts and scn1labΔ44 homozygotes spent less time active compared to wild types (all three experiments p < 0.001 by two-way ANOVA).

is based on ANOVA results computed by Vp_Analyse and found in *190813_both_anova.pdf*.

Figure 6A right (BOX7): `190819_traceByGenotype_BOX6.R` uses *190813_06_deltapxsecsmooth.csv* as input, uploads *190813_06_lbsec.csv* (timepoints of the day/night boundaries) and *080819_06genotype.txt*, and generates *f0_scn1labtrace_07_small.png*.

Figure 6—figure supplement 1 (BOX6): `190819_traceByGenotype_BOX7.R` uses *190813_07_deltapxsecsmooth.csv* as input, uploads *190813_07_lbsec.csv* (timepoints of the day/night boundaries) and *080819_07genotype.txt*, and generates *f0_scn1labtrace_06_small.png*.

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
