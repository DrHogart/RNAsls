#!/usr/bin/perl

# This script searches the hairpins similar to the ILS and GLS - RNA localization
# signals in the I-element and Gurken mRNAs respectively.
# The script is based on the source code from paper of Hamilton et al. 2009
# (http://dx.doi.org/10.1261/rna.1264109).
# The original source code: https://github.com/darogan/RNA2DSearch
#
# The program is distribited under GNU license, v.3:


use strict;
use Getopt::Std;

# fetch options
my %options;
getopts("t:w:i:m:",\%options);

my $T = $options{'t'}; #temp, 25
my $W = $options{'w'}; #window (in nt.), 64
my $max_distance = $options{'m'}; #threshold, 0.11
my $InputFile = $options{'i'}; # Input file

# Canonical RNA loc signals:
my @QuerySeqs = 
   ("aaguaauuuucgugcucucaacaauugucgccgucacagauuguuguucgagccgaaucuuacu",     # GLS
    "ugcacaccucccucgucacucuugauuuuucaagagccuucgaucgaguaggugugca",           # ILS
    "ugcuggcaauccucuccaaaauacucgaaagaguauuucugcggagaguguugccaguac",         # G2LS
    "gugauacgcucgcggauuaacaaggggcucugaggagcuucggauauuccucagcgaucac"         # JLS
   );

my @SS_Seqs   = 
   (".(((((..((((.((((.((((((((((.......))).)))))))...))))))))..)))))",     # GLS
    ".((((((((..((((((.(((((((....))))))).....)).)))).)))))))).",           # ILS
    "(((((((((..((((((.((((((((....))))))))....))))))..))))))))).",         # G2LS
    "(((((.((((.(.(((.......(((((((....)))))))......))).))))))))))"         #JLS
   );


# Parsing of input file
my ($name, $InputFileParsed) = inputFile_parsing($InputFile);

# Run of RNALfold
my $OutputFile  = $InputFileParsed;
$OutputFile =~ s/.no_nl//g;
$OutputFile .= ".output_$T" . "_$W";
print `RNALfold -L $W -T $T < $InputFileParsed > $OutputFile`;

# Parsing of RNALfold results
my @folded_hairpins = RNALfold_parsing($OutputFile, $name);

# measuring of distance between detected haipins and canonical RNA loc signals
my @distances = distance_measure(@folded_hairpins);

# filtering of results with max_distance
my @filtered_distances = distance_filtering($max_distance, @distances);

print STDOUT @filtered_distances;

#### SUBROUTINES HERE #####
#--------------------------

sub inputFile_parsing {

  my($InputFile, $Header, $Length, $WholeInput, $NewInputFile);

  $InputFile = $_[0];
  open(INPUT, "$InputFile") || die "Can't open $InputFile: $!\n";
  while (<INPUT>) { $WholeInput .= $_; }
  close INPUT || die "Can't close $InputFile: $!\n";

  # Deal with header
  $WholeInput =~ m/(>.*)/;
  $Header = $1;
  $Header =~ s/^>//g;
  chomp($Header);

  # Remove Newline characters, for RNALfold input
  $WholeInput =~ s/^\>.*//g;
  $WholeInput =~ s/\n//g;

  $NewInputFile = $InputFile . ".no_nl";

  open (NEWINPUT, ">$NewInputFile") || die "Can't open $NewInputFile: $!\n";
  print NEWINPUT "$WholeInput\n";
  close NEWINPUT || die "Can't close $NewInputFile: $!\n";

  return ($Header, $NewInputFile);
}


sub RNALfold_parsing {

  my $file = @_[0];
  my $name = @_[1];
  my @parsed_results;

  # Deal with individual structure predictions
  open(STRUCTOUT, "$file") || die "Can't open $file: $!\n";
  my @STRUCT_OUT = <STRUCTOUT>;
  close STRUCTOUT || die "Can't close STRUCTOUT: $!\n";

  my $input_sequence = $STRUCT_OUT[@STRUCT_OUT-2];

  for(my $j=0; $j<=$#STRUCT_OUT; $j++) {
    if($STRUCT_OUT[$j] =~ m/^[\.\(\)]/) {
    #Deal With each structure output line from RNALfold
        $STRUCT_OUT[$j] =~ s/\( /\(/g;
        my ($str, $dG, $start) = split(/\s+/, $STRUCT_OUT[$j]);
        # Structure
        $str =~ s/ //g;
        chomp($str);
        # Energy
        $dG =~ s/\(//g;
        $dG =~ s/\)//g;
        chomp($dG);
        # Start
        $start =~ s/ //g;
        chomp($start);
        # Get the sequence for the structure
        my $Seq = substr($input_sequence, $start-1, length($str));
        # Print to big output file
        my $name_subseq = $name.'_'.$j;
        push @parsed_results, join("\t", $name_subseq, $str, $dG, $start, $Seq);
      }
    }

  unlink($file);
  unlink($InputFileParsed);
  return(@parsed_results)
}


sub distance_measure {

  my @Structures = @_;
  my @resulted_distance;

  for(my $i=0; $i<=$#Structures; $i++) {
    my ($ID, $str, $dG, $start, $sequence) = split(/\t/, $Structures[$i]);
    my @DistScore;

    for(my $j=0; $j<=$#QuerySeqs; $j++) {
      # Compare each structure to the query & get a score for it
      my $StructFile = "temp.file";
      open(STRUCT, ">$StructFile") || die "Can't open $StructFile: $!\n";
      print STRUCT "$SS_Seqs[$j]\n$str\n";
      close STRUCT || die "Can't close $StructFile: $!\n";

      my $Dist_Out = `RNAdistance -DfhwcFHWC < $StructFile`;
      $Dist_Out =~ m/f: (.*) h: (.*) w: (.*) c: (.*) F: (.*) H: (.*) W: (.*) C: (.*)/;
      $DistScore[$j] = $1; 
      $DistScore[$j] =~ s/ //g; 
      chomp($DistScore[$j]); 

      unlink($StructFile);
      undef($Dist_Out); 

      #RNAforester similarity
      my $ForesterFile = "temp.file";

      open (FOREST, ">$ForesterFile") || die "Can't open $ForesterFile: $!\n";
      print FOREST ">Query\n$QuerySeqs[$j]\n$SS_Seqs[$j]\n";
      print FOREST ">ID\n$sequence\n$str\n";
      close FOREST || die "Can't close $ForesterFile: $!\n";

      my $Forest_Out_1 = `RNAforester --score --noscale -l < $ForesterFile`;
      chomp($Forest_Out_1); 
      $Forest_Out_1 =~ s/ //g;

      unlink($ForesterFile);
      unlink("rna.ps");
      unlink("x_1.ps");
      unlink("y_1.ps");
 
      $DistScore[$j] .= ',' . sprintf("%.0f", $Forest_Out_1);

      undef($Forest_Out_1); 
    }

    push @resulted_distance, sprintf("REP: ID=%s Length=%-3d Start=%-3d Score=%s-%s-%s-%s dG=%s",
                     $ID, length($sequence), $start, $DistScore[0], $DistScore[1], $DistScore[2], $DistScore[3], $dG);
    push @resulted_distance, "SEQ: $sequence", "SS_: $str";
  }

  undef(@Structures);

  return(@resulted_distance)
}

sub distance_filtering {

  my $threshold = shift @_;
  my @Structures = @_;
  my @filtered_distance;

  # normalization: RNAdistance x/100 and RNAforester 1/x
  for(my $i=0; $i<=$#Structures; $i+=3) {
    my ($rep, $ID, $size, $start, $scores, $dG) = split(/\s+/, $Structures[$i]);
    my $sequence = $Structures[$i+1];
    chomp $sequence;
    $sequence =~ s/SEQ: /Seq=/;
    my $str = $Structures[$i+2];
    $str =~ s/SS_: /Str=/;
    $scores =~ s/Score=//;
    my @all_scores = split('-', $scores);
    my @new_score;
    my $min_value = 1000;
    for (@all_scores) {
      my @s = split(/,/);
      $s[0] = $s[0]/100;
      $s[1] = 2.5/$s[1];
      if ($s[0] < $min_value) {
        $min_value = $s[0];
      }
#      if ($s[1] < $min_value) {
#        $min_value = $s[1]
#      }
      push @new_score, sprintf("[%.2f][%.2f]", $s[0], $s[1]);
    }
#    print STDOUT join(" ", $min_value, $ID, $scores, join(':', @new_score)), "\n";
    # filtering
    if ($min_value <= $threshold) {
      my $scores = join('-', @new_score);
      push @filtered_distance, join("\t", $ID, $start, 'Scores='.$scores, $sequence, $str, $dG),"\n";
    }
  }

  undef(@Structures);
  return(@filtered_distance);

}

