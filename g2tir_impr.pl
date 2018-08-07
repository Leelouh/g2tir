#!/usr/bin/perl
use strict;
use warnings 'all';
use Bio::SeqIO;
use Data::Dumper;
use Bio::DB::Fasta;
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

my $max_length = 4000;
my $min_length = 43;

#my $score = 1;

my $maxScore = 10;

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


my ($motifPrint,$motifDegenerated, $iupac, $fa, $out_fasta, $sum, $fhfasta, $printScore, $TIR1, $TIR2, $score);

my $time = time;

my $output = "g2tir_out.$time";


GetOptions ("in=s"=> \$fa,
            "motif=s"=> \$iupac,
            "out=s"=> \$output,
            "maxL"=> \$max_length,
            "minL"=> \$min_length,
            "pm" => \$motifPrint, #print les TIR
            "degen" => \$motifDegenerated, #dégénère le motif
            "fasta|fa"=> \$out_fasta,
            "config"=> \$sum, #file output config
            "printScore"=> \$printScore,
            "TIR1|R1=s"=> \$TIR1,
            "TIR2|R2=s"=> \$TIR2,
            "score"=> \$score
            );

die "You must specify a fasta file (option --in fasta_file)\n" if ! defined $fa;
die "You must specify a motif in iupac format (option --motif iupac_motif)\n" if ! defined ($iupac || $TIR1);
die "The fasta file <$fa> doesn't exist" if ! -e $fa;
die "Maximum and minimum lengths are wrong, maxL must be higher than minL\n" if $max_length<$min_length;

if ($out_fasta){
  my $output_fasta = "$output.fa";
  open ($fhfasta,'>', $output_fasta) or die "Could not open file '$output_fasta' $!";;
}

open (my $fh, '>', $output);

sub iupacToRegex {
  my ($self) = @_;
  my $re = join '', map $IUPAC{ $_ }, split '', $self;
  return $re;
}

my ($motifDirect, $motifIndirect);

if ($iupac){
  $motifDirect = iupacToRegex($iupac);
}
elsif ($TIR1){
  $motifDirect = iupacToRegex($TIR1);
  $motifIndirect = iupacToRegex($TIR2);
}

#print $motifDirect,"\t",$motifIndirect,"\n";

sub tailleMotif {
  my ($motif) = @_;
  my $cpt_lettre = 0; #initialisation du compteur de la taille du motif
  my @motifD_tab = split //,$motif; #découpage de la string en tableau
  foreach my $l (@motifD_tab){            #on parcourt le tableau
    $cpt_lettre++;                        #à chaque lettre on incrémente le compteur
    if ($l =~ m/\]/){                     #si on rencontre un ] on "décrémente" le compteur
      $cpt_lettre--;
    }
  }
  my $x = "\\[[A-Z]*\\]";                 #on affecte à x un motif en particulier
                                          # le motif en question [[A-Z]*]
  my @c = $motif =~ m/$x/g;         # on récupère dans le tableau chacun des matchs
  my $count = @c;                         # dans $count on retrouve le nombre de case
  foreach my $car (@c){                   # on parcourt chacune des cases
    my @nbc = $car =~ m/[A-Z]/g;          # compte le nombre de lettre à l'intérieur de []
    my $nbc = @nbc;
    $cpt_lettre = $cpt_lettre - $nbc;     #déduction des lettres du compteur de la taille du motif
  }
  return $cpt_lettre;
}

my $motifSizeDirect = tailleMotif($motifDirect);
#my $motifSizeIndirect = tailleMotif($motifIndirect);

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
  my (@patternP, @patternM);
  push (@patternP, $mo);
  my $revmo;
  if ($iupac){
    $revmo = revCompMotif($mo);
  }
  elsif ($TIR1){
    $revmo = $TIR2;
  }
  push (@patternM, $revmo);
  return (\@patternP, \@patternM);
}

#option --degen : si activé, le motif est dégénéré, sinon, non
my ($tabP, $tabN);

if ($motifDegenerated){
  if ($iupac){
    my $revMotif = revCompMotif($motifDirect);
    $tabP = motifDegenere($motifDirect);
    $tabN = motifDegenere($revMotif);
  }
   elsif ($TIR1){
#     ($tabP, $tabN) = motifDegenere($motifDirect,$motifIndirect);
    $tabP = motifDegenere($motifDirect);
    $tabN = motifDegenere($motifIndirect);
   }
}
else {
  if ($iupac){
    ($tabP, $tabN) = notDegenerated($motifDirect);
  }
  elsif ($TIR1){
    ($tabP, $tabN) = notDegenerated($motifDirect);
  }
}
#print $tabP,"\t",$tabN,"\n";
my $try = findPairsMotif($tabP, $tabN);

sub findPairsMotif {
  my ($tabP, $tabN)= @_;
#  print @$tabP,"\t",@$tabN,"\n";
  my $seqio_obj = Bio::SeqIO->new( -file => "$fa",
                      -format => "fasta" );

  my $listP = join ('|', @$tabP);
  my $listN = join ('|', @$tabN);

  while ( my $seq_obj = $seqio_obj->next_seq ) {
    my $seq = $seq_obj-> seq;
    my $chr = $seq_obj->display_id;
    my @endPos;
    my $maxSkyline;
    my $minSkyline;
    my $startSubstr;
    my $length_extension;

    while ($seq =~ m/$listP/gio){
      $chr =~ s/(chr* ).*/$1/;
      my $pPos = pos($seq);
      my $pStart = $pPos - $motifSizeDirect;
      my $pEnd = $pPos;
      if (@endPos){
        my $index = 0;
        for my $value (@endPos){
          if ($value < ($pEnd + $min_length + $motifSizeDirect)){
            $index++;
          }
          else {
	           last
          }
        }
        @endPos = $index < scalar @endPos ? @endPos[$index..$#endPos] : ();
      }
      $minSkyline = $pEnd + $min_length;
      my $R1 = substr($seq, $pStart, $motifSizeDirect);
      # print "startSub : $startSubstr\n";
      # print "R1 : $R1\n";
      # print "pStart : $pStart\n";
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
      # print "sous chaine : $substr\n";
      # print "LISTN : $listN\n";
      while ($substr =~ m/$listN/gio){
      	my $nPos = pos($substr);
      	my $nEnd = $nPos + $startSubstr;
      	push @endPos, $nEnd;
      }

      my $motifSizeIndirect;
      if ($iupac){
        $motifSizeIndirect = $motifSizeDirect;
      }
      elsif ($TIR1){
        $motifSizeIndirect = tailleMotif($motifIndirect);
      }
      foreach my $nEndPos (@endPos){
      	my $nStart = $nEndPos - $motifSizeIndirect;
      	my $R2 = substr($seq, $nStart, $motifSizeIndirect);
      	my $R2_uc = revCompMotif(uc($R2));

      	my $dist = distance($R1_uc, $R2_uc);
        #my @out;

        my $scoreComp;

        if ($score){
          $scoreComp = $score;
        }
        else {
          $scoreComp = $maxScore;
          if ($R1){
            $scoreComp = 200;
          }
        }

        if ($iupac){
            if ($dist <= $scoreComp){
    #    my @out = ($chr, $pStart, $nEndPos, $nEndPos- $pStart, $dist);
    #          $pStart=$pStart+1; #mise à niveau pour standart extraction fasta
              my @out = ($chr, $pStart, $nEndPos);
              push @out, ($R1, $R2) if $motifPrint;
              push @out, ($dist) if $printScore;
              print $fh join ("\t", @out), "\n";
              if ($out_fasta){
    #            $pStart =$pStart+1;
                my $PBLE = $seq_obj->subseq($pStart+1, $nEndPos);
                print $fhfasta ">$chr:$pStart-$nEndPos\n$PBLE\n";
              }
            }
          }

        elsif ($R1){
          if ($dist <= $scoreComp){
            my @out = ($chr, $pStart, $nEndPos);
            push @out, ($R1, $R2) if $motifPrint;
            push @out, ($dist) if $printScore;
            print $fh join ("\t", @out), "\n";
            if ($out_fasta){
  #            $pStart =$pStart+1;
              my $PBLE = $seq_obj->subseq($pStart+1, $nEndPos);
              print $fhfasta ">$chr:$pStart-$nEndPos\n$PBLE\n";
            }
          }
        }

      }
    }
  }
}

if ($sum){
  my $file_sum = "$output.config";
  open (my $fhsum, '>', $file_sum);
  print $fhsum "motif : $iupac\nfasta : $fa\n" if ($iupac);
  print $fhsum "motif was degenerated\n" if ($motifDegenerated);
  print $fhsum "fasta file was generated\n" if ($out_fasta);
}
