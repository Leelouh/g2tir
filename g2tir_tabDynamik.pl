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
# perl test.pl motif genome
# ex : perl test.pl TTAAC[CAT][CT].[TCA] /home/lhelou/data/human_data/grch38-p10.fa

#my $output = "/home/lhelou/Documents/g2tir/test/chr22_xsCol.out";

my $output = "chr2_dynamik_test.out";


#my $output = "/home/lhelou/scripts/perl/g2TIR/new/test/chr22.out";
open (my $fh, '>', $output);


#my $fa = "/home/lhelou/Documents/g2tir/chr1_ech_test.fa";
#my $fa = "/home/lhelou/scripts/perl/g2TIR/new/test/chr22.fa";

my $fa = "/home/lhelou/Data/chr22.fa";
$fa = '/data/dbseq/Homo_sapiens/hg38/chr/chr22.fa';
#$fa = 'chr1_2_test.fa';

#my $fa = "/home/lhelou/data/human_data/grch38-p10.fa";

my $max_length = 4000;
my $min_length = 43;

##VRAI
#my $motifDirect = "TTAAC[CAT][CT].[TCA]";
#my $motifRegex = "TTAACHYNH";
#my $motifDirect = iupac($motifRegex);
#print $motifDirect;

#my $motifDirect = "TTAACTTTT";

my ($motifPrint,$motifDegenerated);

GetOptions ("motif" => \$motifPrint, #print les TIR
            "degen" => \$motifDegenerated #dégénère le motif
            );

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



sub iupacToRegex {
  my ($self) = @_;
  #print "self : $self\n";
  my $re = join '', map $IUPAC{ $_ }, split '', $self;
  #print $re,"\n";
  return $re;
}

# print $motifDirect,"\n";

my $iupac = "TTAACHYNH";
my $motifDirect = iupacToRegex($iupac);
#my $motifDirect = "TTAAC[CAT][CT].[TCA]";


#my $motifDirect = $ARGV[0];
#my $fa = $ARGV[1];

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
  #print $cpt_lettre;
  return $cpt_lettre;
}

my $motifSize = tailleMotif($motifDirect);
#my $lengthDiv = tailleMotif($motifDirect)-1;

sub revCompMotif { # fonction permettant de générer le reverse comp du motif étudié
#print "revcom\n";
  my ($motif) = @_;# possible de reverse .{2} mais pas possible de reverse [AC]{2} => need amélioration
  my $rev = reverse($motif);
  $rev =~ tr/ACGT[]/TGCA][/;
  $rev =~ s/\}(.*)\{\./\.{$1\}/;
  return $rev;
}


sub motifDegenere {
  #print "motifDegenere\n";
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

    if ( ($ss ne "]" ) && ($ss ne "[")){
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

  for my $pos (@patternP){
    my $motI=revCompMotif($pos);
    push (@patternM, $motI);
  }

  # foreach my $ligne (@patternP){
  #   print $ligne,"\n";
  # }
  # foreach my $ll (@patternM){
  #   print $ll,"\n";
  # }
  return (\@patternP, \@patternM);
}

sub notDegenerated {
  my ($mo) = @_;
  my (@patternP, @patternM);
  push (@patternP, $mo);

  my $revmo = revCompMotif($mo);
  push (@patternM, $revmo);

  return (\@patternP, \@patternM)
}

#option --degen : si activé, le motif est dégénéré, sinon, non
my ($tabP, $tabN);
if ($motifDegenerated){
  ($tabP, $tabN) = motifDegenere($motifDirect);
}
else {
  ($tabP, $tabN) = notDegenerated($motifDirect);
}

my $try = findPairsMotif($tabP, $tabN);


sub findPairsMotif {

  my ($tabP, $tabN)= @_;
  # foreach my $ligne (@$tabN){
  #   print $ligne,"\n";
  # }
  my $db = Bio::DB::Fasta -> new($fa);
  my $seqio_obj = Bio::SeqIO->new(-file => "$fa", #on recupere la sequence depuis notre fichier (arg 0)
                      -format => "fasta" );

  my $listP = join ('|', @$tabP);
  my $listN = join ('|', @$tabN);

  while ( my $seq_obj = $seqio_obj->next_seq ) {
    #my %h;
    my $seq = $seq_obj-> seq;
    my $chr = $seq_obj->display_id;

    my @startPos;
    my @endPos;
    my $maxSkyline;
    my $minSkyline;
    my $startSubstr;
    my $length_extension;

    while ($seq =~ m/$listP/gio){
      $chr =~ s/(chr* ).*/$1/;
      my $pPos=pos($seq);
      my $pStart = $pPos - $motifSize;
      my $pEnd = $pPos;

      #warn Dumper \@endPos;
      #warn sprintf "term: %i\n", $pEnd + $min_length;
      if (@endPos){
        my $index = 0;
        for my $value (@endPos){
          if ($value < ($pEnd + $min_length)){
	    $index++;
          } else {
	    last
	  }
        }
	@endPos = $index < scalar @endPos ? @endPos[$index..$#endPos] : ();
      }
      #warn Dumper \@endPos;
      #warn sprintf "$chr term: %i -- %i elem -- range=%i (%i:%i)\n", $pEnd + $min_length, scalar @endPos, 
      #scalar @endPos ? ($endPos[$#endPos]- $endPos[0]+1,$endPos[0], $endPos[$#endPos]) : (0, 0, 0);
      $minSkyline = $pEnd + $min_length;
      $maxSkyline = $minSkyline + $max_length;

      my $R1 = substr($seq,$pStart,$motifSize);
      my $R1_uc = uc($R1);
      if (@endPos){
	$startSubstr = $endPos[-1];
	$length_extension = $max_length-($startSubstr-$pEnd - $min_length);
      }
      else {
	$startSubstr = $minSkyline;
	$length_extension = $max_length;
      }
      #warn sprintf ("# substr %i %ib\n", $startSubstr, $length_extension);
      my $substr = substr($seq, $startSubstr,$length_extension);

      while ($substr =~m/$listN/gio){
	my $nPos = pos($substr);
	my $nEnd = $nPos+$startSubstr;
	push @endPos, $nEnd;
      }
      foreach my $nEndPos (@endPos){
	my $nStart=$nEndPos-$motifSize;
	my $R2=substr($seq,$nStart, $motifSize);
	my $R2_uc=revCompMotif(uc($R2));

	my $dist = distance($R1_uc,$R2_uc);
	if ( $dist <= 1){
	  my @out = ($chr, $pStart, $nEndPos, $nEndPos- $pStart, $dist);
	  push @out,  ($R1, $R2) if $motifPrint;
	  print $fh join ("\t", @out), "\n";
	}
      }
    }
  }
}

          __END__
          while ($substr  =~m/$listN/gio){

            $nPos = pos($substr);
            # my $nStart=$nPos-$motifSize+$minSkyline;
            # my $nEnd = $nPos+$minSkyline;
            #$nStart=$nPos-$motifSize+$length_extension;
            my $nEnd = $nPos+$length_extension;
#            print "\n$pStart\t$nStart\t$startSubstr\t$length_extension=====\n";

            #my $R2=substr($seq,$nStart, $motifSize);

            # my $R1_uc = uc($R1);
            # my $R2_uc=revCompMotif(uc($R2));

            push (@endPos,$nEnd);
            #my $PBLE = substr($seq,$pStart,($nEnd-$pStart));
            #my $distances = distance($truc,$revsubC);
            #print "$chr\t$start\t$endNe\t$truc\t$revsub\t$distances\n";

          }
          foreach my $endMinus (@endPos){
            $nStart=$endMinus - $motifSize;
            # my $R2=substr($seq,$nStart, $motifSize);
            # my $R1_uc = uc($R1);
            # my $R2_uc=revCompMotif(uc($R2));
            #
            # if ( distance($R1_uc,$R2_uc) <= 1){
            #   my $debT = $pStart+1; #sinon samtools faidx commence une base plus tôt (décalage)
            #         #  print "R1 : $motif1\nR2 : $motif2Temp\nREV: $motif2\n\n";
            #         #  print $a,"\n",$b,"\n\n";
            #   print "$chr\t$debT\t$endMinus\t",distance($R1_uc,$R2_uc);
            #   if ($motifPrint){
            #     print "\t$R1\t$R2";
            #   }
            #   print "\n";
            # }
          }
            # if ( distance($R1_uc,$R2_uc) <= 1){
            #     my $debT = $pStart+1; #sinon samtools faidx commence une base plus tôt (décalage)
            #       #  print "R1 : $motif1\nR2 : $motif2Temp\nREV: $motif2\n\n";
            #       #  print $a,"\n",$b,"\n\n";
            #       push (@endPos, $nEnd);
            # print "$chr\t$debT\t$nEnd\t",distance($R1_uc,$R2_uc);
            # if ($motifPrint){
            #   print "\t$R1\t$R2";
            # }
            # print "\n";
            # }
         #}
        }
      }
  }
}
