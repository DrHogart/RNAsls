#!/usr/bin/perl

# This script searches the hairpins similar to the ILS and GLS - RNA localization
# signals in the I-element and Gurken mRNAs respectively.
# The script is based on the source code from paper of Hamilton et al. 2009
# (http://dx.doi.org/10.1261/rna.1264109).
# The original source code: https://github.com/darogan/RNA2DSearch
#
# The program is distribited under GNU license, v.3.


use strict;
use Getopt::Std;
use Bio::SeqIO;

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


my @resulted_table;

# Loading of input file
my $seq_in = Bio::SeqIO->new(-file => "<$InputFile", 
                             -format => 'fasta');

# Deal with individual sequences: find the candidate hairpins
while (my $seq = $seq_in->next_seq) {
  
  # the name of the current sequence
  my $name = $seq->primary_id();
  
  # call RNALfold
  my $RNALfold_file = run_RNALfold($seq);

  # Parsing of RNALfold results
  my @folded_hairpins = RNALfold_parsing($RNALfold_file, $name); 

  # measuring of distance between detected hairpins and canonical RNA loc signals
  my @distances = distance_measure(@folded_hairpins);

  # discarding candidate hairpins having too large distance (> $max_value)
  my @filtered_distances = distance_filtering($max_distance, @distances);

  # calculate randfold p-value
  my @add_randfold = compute_randfold_pvalue(\@filtered_distances);
  print @add_randfold;

  # collect the results
  push @resulted_table, @filtered_distances;

}

#print STDOUT @resulted_table;


#### SUBROUTINES HERE #####
#--------------------------

sub run_RNALfold {

  my $seq = @_[0];

  # rand number for temp files
  my $rand = int(rand(999));
  my $temp_fasta = 'temp.'.$rand.'.fasta';
  my $temp_str = 'temp.'.$rand.'.str';

  my $seq_out = Bio::SeqIO->new(-file   => ">$temp_fasta",
                                -format => 'fasta');
  $seq_out->write_seq($seq);

  # Run of RNALfold
  print `RNALfold -L $W -T $T < $temp_fasta > $temp_str`;

  unlink($temp_fasta);
  return($temp_str);

}


sub RNALfold_parsing {

  my $file = @_[0];
  my $name = @_[1];
  my @parsed_results;

  open(STRUCTOUT, "$file") || die "Can't open $file: $!\n";
  my @strs = <STRUCTOUT>;
  close STRUCTOUT || die "Can't close STRUCTOUT: $!\n";

  # fetch the whole sequence
  my $input_sequence = $strs[@strs-2];

  # Deal with individual structure predictions
  for(my $j=0; $j<=$#strs-1; $j++) {
    if($strs[$j] =~ m/^[\.\(\)]/) {
        #Deal With each structure output line from RNALfold
        chomp $strs[$j];
        $strs[$j] =~ s/\( /\(/g;
        my ($str, $dG, $start) = split(/\s+/, $strs[$j]);
        $str =~ s/ //g;
        $start =~ s/ //g;
        $dG =~ s/[\(\)]//g;
        # Get the sequence for the structure and collect the results
        my $Seq = substr($input_sequence, $start-1, length($str));
        push @parsed_results, join("\t", $name.'_'.$j, $str, $dG, $start, $Seq);
      }
    }

  unlink($file);
  return(@parsed_results)
}


sub distance_measure {

  my @strs = @_;
  my @resulted_distance;

  for(my $i=0; $i<=$#strs; $i++) {
    my ($ID, $str, $dG, $start, $sequence) = split(/\t/, $strs[$i]);
    my @DistScore;

    for(my $j=0; $j<=$#QuerySeqs; $j++) {
      # Compare each structure to the query & get a score for it
      my $rand = int(rand(999));
      my $StructFile = "temp.".$rand.".dat";
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
      my $ForesterFile = "temp.".$rand.".dat";

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

      unlink($StructFile);
      undef($Forest_Out_1); 
    }

    push @resulted_distance, sprintf("REP: ID=%s Length=%-3d Start=%-3d Score=%s-%s-%s-%s dG=%s",
                     $ID, length($sequence), $start, $DistScore[0], $DistScore[1], $DistScore[2], $DistScore[3], $dG);
    push @resulted_distance, "SEQ: $sequence", "SS_: $str";
  }

  undef(@strs);
  return(@resulted_distance)
}

sub distance_filtering {

  my $threshold = shift @_;
  my @strs = @_;
  my @filtered;

  # normalization: RNAdistance x/100 and RNAforester 1/x
  for(my $i=0; $i<=$#strs; $i+=3) {
    my ($rep, $ID, $size, $start, $scores, $dG) = split(/\s+/, $strs[$i]);
    my $sequence = $strs[$i+1];
    chomp $sequence;
    $sequence =~ s/SEQ: /Seq=/;
    my $str = $strs[$i+2];
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
      push @new_score, sprintf("[%.2f][%.2f]", $s[0], $s[1]);
    }

    if ($min_value <= $threshold) {
      my $scores = join('-', @new_score);
      push @filtered, join("\t", $ID, $start, 'Scores='.$scores, $dG, $sequence, $str)."\n";
    }
  }

  undef(@strs);
  return(@filtered);

}

sub compute_randfold_pvalue {

  my $ref_rnals = $_[0];
  my @added_pvalues;

  for my $s (@$ref_rnals) {
    chomp $s;

    my ($id, $start, $scores, $dG, $sequence, $str) = split(/\t/, $s);
    $sequence =~ s/Seq=//;
    
    my $rand = int(rand(999));
    my $fasta_file = 'temp.'.$rand.'.fasta';
    open (my $fh, ">$fasta_file");
    print $fh '>'.$id."\n".$sequence."\n";
    close $fh;

    my $cmd_randfold = `randfold -d $fasta_file 999`;
    chop $cmd_randfold;
    my $pvalue = (split(/\s+/, $cmd_randfold))[2];
    push @added_pvalues, join("\t", $id, $start, $scores, $dG, 'P='.$pvalue, 'Seq='.$sequence, $str)."\n";

    unlink($fasta_file);
  }
  return @added_pvalues;
}

