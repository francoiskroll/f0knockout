# *trpa1b* mustard oil assay

Refers to Figure 5A.

___

Please cite the preprint if you use some of the data or code! <br />
https://doi.org/10.1101/2020.06.04.133462

___

*070220_07_genotypeMap.xlsx*: layout of the plate recorded in the ZebraBox. <br />
Converted into genotype file *070220_07genotype.txt* with `genotypeGenerator.R`.

At the end of the experiment, ZebraLab generates a raw file. <br />
In ZebraLab software: it is exported into xls files (each 50,000 rows), with: <br />
Raw data... > Export > ...

The xls files are in *200211_06_07_mustard_rawoutput.zip*.

The genotype file + xls files are input into MATLAB script `VpExtract_ErrorHandling.m` (in folder **behaviour_matlabscripts**/**Vp_Extract**), originally forked from https://github.com/ghoshm/Structure_Paper <br />
(Essentially using the first part which handles the ZebraBox tracking glitches and appends the xls files).

`VpExtract_ErrorHandling.m` generates *070220_07_deltapxsqsec.csv* (included).

*070220_07_deltapxsqsec.csv* is then used as input to R script `traceByGenotypeSMALL.R`. <br />
`traceByGenotypeSMALL.R` generates the pre/post activity traces of Figure 4A left.

*070220_07_deltapxsqsec.csv* is then used again as input to R script `trpa1b_lineplots.R`. <br />

`trpa1b_lineplots.R`:
* generates the line plot in Figure 4A right
* summarises the results for each fish in delta px (total px post - total px pre), and exports the data to *dpxresults.csv* (included).
* does a t-test on the delta px values.

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
