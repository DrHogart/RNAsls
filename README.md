
# RNAsls: the search of RNA localization signals


The RNAsls searches the hairpins similar to the ILS and GLS - RNA localization signals in the I-element and Gurken mRNAs respectively. These signals are necessary for the proper transport of RNA along microtubes during Drosophila oogenesis.
The code of RNAsls is just a scaled-down version of RNA2DSearch program published by Hamilton et al. 2009 (http://dx.doi.org/10.1261/rna.1264109).
The original source code of RNA2DSearch: https://github.com/darogan/RNA2DSearch.

Usage
=====
* Prerequisites: [BioPerl](http://www.bioperl.org/wiki/Main_Page), [RNA Vienna package](http://rna.tbi.univie.ac.at), [RNAforester](http://bibiserv.techfak.uni-bielefeld.de/rnaforester/), [randfold](http://bioinformatics.psb.ugent.be/supplementary_data/erbon/nov2003).

* Example of usage:
```html
perl RNAsls.pl -i InputFile.fasta -t 25 -w 64 -m 0.11 > OutputFile.txt
```

* Options
    - `-i` Input file with sequence(-s) in fasta format
    - `-t` Temperature for RNA folding, in Celsius scale (e.g. 25)
    - `-w` Maximum allowed length of the candidate hairpins, in nt. (e.g. 64)
    - `-m` Maximum allowed distance between candidate hairpin and the canonical RNA localization signals (e.g. 0.11). 'Distance' is a normalized [RNAdistance](http://rna.tbi.univie.ac.at/cgi-bin/RNAfold.cgi) score.

* Output: The table with subsequences that have secondary structure similar to the canonical RNA localization signals. The level of similarity (distance) is presented in the `Scores` field as the set of score pairs. The first score in each pair is normalized [RNAdistance](http://rna.tbi.univie.ac.at/cgi-bin/RNAfold.cgi) score, while the second one is normalized and inversed [RNAforester](http://bibiserv.techfak.uni-bielefeld.de/rnaforester/) score. Since the comparison is performed with four known RNA localization signals, including GLS, ILS, G2LS, and JLS, the `Scores` field contains four pairs of scores. The output table also includes the `P`-value of [randfold](http://bioinformatics.psb.ugent.be/supplementary_data/erbon/nov2003) test of hairpin stability.

License
=======

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
