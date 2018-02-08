Versions
========

* 2.1 (2018-02-08)
  	* Bug fixes and performance improvements.
	* New program: iu-gen-matching-fastq-files. Helps avoiding getting cancer when people release R1/R2 files that have different number of reads and non-matching pairs.

* 1.4.8 (2016-06-10)
	* iu-merge-pairs bug fixed: when `--marker-gene-stringent` is used, the program no longer assumes that the user wants to retain only the overlapping parts of reads (for which we have a separate flag: `--retain-only-overlap`).

* 1.4.7 (2016-05-11)
	* iu-trim-V6-primers no longer trims information-rich deflines

* 1.4.6 (2016-05-10)
	* Realizing that there was this CHANGELOG.
	* Improvements over pretty much everything.

* 1.4.1 (2014-05-27)
    * Better error reporting.

* 1.3 (2015-02-25)
    * --ignore-defline parameter has been added to Minoche and Bokulich QC.

* 1.2 (2015-02-04)
    * Partial support for very old versions of CASAVA

* 1.1 (2015-01-27)
    * Scripts are renamed (now they all have 'iu-' prefix).
    * Minor enhancements.

* 1.0 (2015-01-24)
    * First stable release.
