#!/usr/bin/perl
use strict;
use warnings 'all';
#use Bio::SeqIO;
#use Data::Dumper;
#use Bio::DB::Fasta;
use Bio::Perl;
#use Text::Levenshtein::XS qw/distance/;
use Text::Levenshtein::XS qw/distance/;
#use Text::Levenshtein::Damerau::XS qw/xs_edistance/;
#use Text::Levenshtein qw(distance);

use Getopt::Long qw/:config bundling auto_abbrev permute/;

# lancement/usage :
# perl g2tir_tabDynamik.pl --pm --degen --motif motif_iupac --fa fasta_file
# ex : time perl g2tir_tabDynamik.pl --pm --degen --motif TTAACHYNH --fa /home/lhelou/Data/chr22.fa

# test : time perl g2tir_tabDynamik.pl --degen --pm --motif TTAAYYYNNNTCATTCCT --in /home/lhelou/data/human_data/grch38-p10.fa --config --fasta --out TTAAYYYNNNTCATTCCT_12_06_18

#my $fa = $ARGV[0];

my $TSD = "TTAA";

my %IUPAC = (
 A => 'A',
 C => 'C',
 G => 'G',
 T => 'T',
 R => '[AG]',
 Y => '[CT]',
 M => '[AC]',
 K => '[GT]',
 W => '[AT]',
 S => '[GC]',
 B => '[CGT]',
 D => '[AGT]',
 H => '[ACT]',
 V => '[ACG]',
 N => '[ACGT]',
);

my $motifSize = 4;

my ($motifPrint,$motifDegenerated, $iupac, $fa, $out_fasta, $sum, $fhfasta);

my $time = time;

my $output = "TSD_search_out.$time";

GetOptions ("in=s"=> \$fa,
            "motif=s"=> \$iupac,
            "out=s"=> \$output,
            "pm" => \$motifPrint, #print les TIR
            "degen" => \$motifDegenerated, #dégénère le motif
            "fasta|fa"=> \$out_fasta
            );

open my $fh, '<', $fa or die;

#print $fa,"\n";

sub extractRegion {
  my ($self) = @_;# split /"\t"/, @_;
  while ($self){
  my ($chr, $start, $end) = split /"\t"/, $self;
  print $chr;
}
}

extractRegion(\$fh);



__END__

sub revCompMotif { # fonction permettant de générer le reverse comp du motif étudié
  my ($motif) = @_;
  my $rev = reverse($motif);
  $rev =~ tr/ACGT[]/TGCA][/;
  $rev =~ s/\}(.*)\{\./\.{$1\}/;
  return $rev;
}

sub motifDegenere {
  my $cpt = 0;
  my ($pattern)=@_;
  my %patterns;
  my %pattAll;
  my (@patternP, @patternM);
  for (my $i=0; $i<length($pattern); $i++){
    my $id = $cpt;
    $cpt++;
    my $copy = $pattern;
    my $ss = substr ($pattern,$i,1);

    if ( ($ss ne "]") && ($ss ne "[") ){
      substr ($copy,$i,1) = ".";
      push (@patternP, $copy);
    }
    else {
      $i++;
      $ss = substr($pattern, $i, 1);
      my @t;
      while($ss ne "]") {
        push(@t, $i);
        $i++;
        $ss = substr($pattern, $i, 1);
     }
    substr ($copy,$t[0],scalar(@t)) = ".";
    push (@patternP, $copy);
    }
  }
  # for my $pos (@patternP){
  #   my $motI=revCompMotif($pos);
  #   push (@patternM, $motI);
  # }
  return (\@patternP);#, \@patternM);
}


sub notDegenerated {
  my ($mo) = @_;
  my (@patternP);
  push (@patternP, $mo);
  return (\@patternP);
}

my ($tab);

if ($motifDegenerated){
  ($tab) = motifDegenere($TSD);
}
else {
  ($tab) = notDegenerated($TSD);
}

#my $motifDegenerated = motifDegenere($TSD);

my $try = findTSD($tab);

sub findTSD {
  my ($tabP)= @_;
  my $seqio_obj = Bio::SeqIO->new( -file => "$fa",
                      -format => "fasta" );

  my $listP = join ('|', @$tabP);

  while ( my $seq_obj = $seqio_obj->next_seq ) {
    my $seq = $seq_obj-> seq;
    my $id = $seq_obj->display_id;
    my @startPos;
    my @endPos;
    my $maxSkyline;
    my $minSkyline;
    my $startSubstr;
    my $length_extension;

    while ($seq =~ m/$listP/gio){
      my $chr = $id;
      $chr =~ s/(chr.*):.*/$1/;
#      my $id =~ s/chr.*//;
      my $pPos = pos($seq);
      my $pStart = $pPos - $motifSize;
      my $pEnd = $pPos;
      my $R1 = substr($seq, $pStart, $motifSize);
      print "séquence : $seq\t$R1\t$pStart\t$motifSize\n";
      ####DEBUGprint "new line : \nid : $id\nchr : $chr\n";
    }
  }
}

      __END__
      if (@endPos){
        my $index = 0;
        for my $value (@endPos){
          if ($value < ($pEnd + $min_length + $motifSize)){
            $index++;
          }
          else {
	           last
          }
        }
        @endPos = $index < scalar @endPos ? @endPos[$index..$#endPos] : ();
      }
      $minSkyline = $pEnd + $min_length;
      my $R1 = substr($seq, $pStart, $motifSize);
      my $R1_uc = uc($R1);
      if (@endPos){
      	$startSubstr = $endPos[-1];
      	$length_extension = $max_length - ($startSubstr - $pEnd - $min_length);
      }
      else {
      	$startSubstr = $minSkyline;
      	$length_extension = $max_length - $min_length;
      }
      my $substr = substr($seq, $startSubstr, $length_extension);
      while ($substr =~ m/$listN/gio){
      	my $nPos = pos($substr);
      	my $nEnd = $nPos + $startSubstr;
      	push @endPos, $nEnd;
      }
      foreach my $nEndPos (@endPos){
      	my $nStart = $nEndPos - $motifSize;
      	my $R2 = substr($seq, $nStart, $motifSize);
      	my $R2_uc = revCompMotif(uc($R2));

      	my $dist = distance($R1_uc, $R2_uc);
      	if ($dist <= 1){
      #    my @out = ($chr, $pStart, $nEndPos, $nEndPos- $pStart, $dist);
          $pStart=$pStart+1;
          my @out = ($chr, $pStart, $nEndPos);
          push @out, ($R1, $R2) if $motifPrint;
          print $fh join ("\t", @out), "\n";
          if ($out_fasta){
            my $PBLE = $seq_obj->subseq($pStart, $nEndPos);
            print $fhfasta ">$chr:$pStart-$nEndPos\n$PBLE\n";
          }
        }
      }
    }
  }
}
