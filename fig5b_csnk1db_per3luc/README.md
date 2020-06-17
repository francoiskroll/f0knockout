# *csnk1db* (*per3:luciferase* assay)

Refers to Figure 5B.

___

Please cite the preprint if you use some of the data or code! <br />
https://doi.org/10.1101/2020.06.04.133462

___

*1689_csv.zip*: csv files obtained from TopCount plate reader (one file per timepoint, each 9.92 min).

Command line script `lucMerge.command` collates the csv files into *220719_all.csv* (included).

*220719genotypeMap.xlsx*: plate layout. Used as input to `genotypeGenerator.R` to generate *220719genotype.txt* (included).

Data file *220719_all.csv* and metadata file *220719genotype.txt* are used as input to script `per3Export.R`. It generates files:

* *220719_per3ts.csv*: the timestamp of each timepoint (not used further)
* *220719_per3cps.csv*: CPS values, with columns names as well/fish number (not used further)
* *220719_per3biodare.csv*: CPS values (same as previous), with columns names as condition; allows BioDare2 to group the data in conditions by reading the header row.

Go to https://biodare2.ed.ac.uk/

Time series file: Coma-separated <br />
Load `220719_per3biodare.csv` (drag and drop) >> Upload selected <br />
Data in columns // Import labels from row // Background noise in columns <br />
Time unit: time in hours <br />
Pick columns as asked <br />
Type to select background label...: "blank" <br />
Import data <br />
Note: it mentions number of timeseries imported; should be number of fish (blank are not counted)

Refer to Materials and Methods for period analysis with BioDare2.

From BioDare2 is exported period analysis results: *dmso_periods.csv* and *pf67_periods.csv*, which are collated into *220719_periodBoth.csv*. <br />
*220719_periodBoth.csv* is used as input to script `per3Period.R`, which generates the period scatterplot in Figure 5B below and does the statistical tests.

From BioDare2 (download the current view) is exported the detrendred/normalised timeseries: *ampbaselinedtr_normmean.csv* (included). <br />
*ampbaselinedtr_normmean.csv* is used as input to script `per3Trace.R`, which generates the timeseries plot in Figure 5B above.

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
