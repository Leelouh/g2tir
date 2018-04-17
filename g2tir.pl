#!/usr/bin/perl
use strict;
use warnings 'all';
use Bio::SeqIO;
use Data::Dumper;
use Bio::DB::Fasta;
use Bio::Perl;
#use Text::Levenshtein::XS qw/distance/;
use Text::LevenshteinXS qw/distance/;
#use Text::Levenshtein::Damerau::XS qw/xs_edistance/;
#use Text::Levenshtein qw(distance);
#use Devel::NYTProf;

use Getopt::Long qw/:config bundling auto_abbrev permute/;

# lancement/usage :
# perl test.pl motif genome
# ex : perl test.pl TTAAC[CAT][CT].[TCA] /home/lhelou/data/human_data/grch38-p10.fa

#my $output = "/home/lhelou/Documents/g2tir/test/chr22_xsCol.out";
#my $output = "chr1_2_test.out";
my $output = "chr22_gtir_test_16_4_18_b.out";
#my $output = "/home/lhelou/Documents/g2tir/test/chr1_ech_degen.out";


#my $output = "/home/lhelou/scripts/perl/g2TIR/new/test/chr22.out";
open (my $fh, '>', $output);


#my $fa = "/home/lhelou/Documents/g2tir/chr1_ech_test.fa";
#my $fa = "/home/lhelou/scripts/perl/g2TIR/new/test/chr22.fa";

#my $fa = "chr1_2_test.fa";
my $fa = "/home/lhelou/Documents/g2tir/echPos.fa";
#$fa = '/data/dbseq/Homo_sapiens/hg38/chr/chr22.fa';
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
  print "self : $self\n";
  my $re = join '', map $IUPAC{ $_ }, split '', $self;
  #print $re,"\n";
  return $re;
}

# print $motifDirect,"\n";

my $iupac = "TTAACHYNH";
my $motifDirect = iupacToRegex($iupac);
#my $motifDirect = "TTAAC[CAT][CT].[TCA]";

#print $ttt;

# GetOptions ("help|?|h" => \$help,
#             "debug+" => \$debug,
#             "only_debug" => \$only_debug,
#             #"mite_name" => \$mite_name,
#             "base_name=s" => \$base_name,
#             "dir=s" => \$out_dir,
#             "bed" => \$bed_format,
#             "rm" => \$rm_input,
#             "tag=s" => \$tag,
#             "annot=s" => \$annot
#             #"trsp=s" => \$str_trsp #exemple : HSMAR1,1287,36
#             );


#my $fa = '/home/lhelou/data/human_data/grch38-p10.fa';

#my $motifDirect = $ARGV[0];
#my $fa = $ARGV[1];

sub tailleMotif {
  my ($motif) = @_;
  my $cpt_lettre = 0; #initialisation du compteur de la taille du motif
  my @motifD_tab = split //,$motif; #découpage de lacoucou string en tableau
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

#my $t = revCompMotif("AATT");
#print $t;
#__END__

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

# sub motifSingle {
#   #print "motifDegenere\n";

#   my ($motif)=@_;
#   my $cpt = 0;
#   my @listMotifs;
# #  print length($motif);
#   for (my $i=0; $i<length($motif); $i++){
#     my $id = $cpt;
#     $cpt++;
#     my $copy = $motif;
#       substr ($copy,$i,1) = ".";
#       push (@listMotifs, $copy);
#   }

#   # foreach my $t(@listMotifs){
#   #   print $t,"\n";
#   # }
#   return (\@listMotifs);
# }


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
    my %h;
    my $seq = $seq_obj-> seq;
    my $chr = $seq_obj->display_id;
    my ($start, $end, $subSeq);



    while ($seq =~ m/$listP/gio){
      $chr =~ s/(chr* ).*/$1/;
        #print $chr;
      my $pos=pos($seq); #pos en partant de la fin de notre séquence


      my $start = $pos - $motifSize; #$lengthDiv;#- 8;
      my $end = $pos;
      my $deb = $end+$min_length;

      my $R1 = substr($seq,$start,$motifSize);
      my $R1_uc = uc($R1);

      my $substr = substr($seq, $deb,$max_length);

      while ($substr  =~m/$listN/gio){

	my $posNe = pos($substr);
	my $start2=$posNe-$motifSize+$deb;
	my $endNe = $posNe+$deb;

	my $R2=substr($seq,$start2, $motifSize);

	my $R2_uc=revCompMotif(uc($R2));

	my $PBLE = substr($seq,$start,($endNe-$start));

	if ( distance($R1_uc,$R2_uc) <= 1){
	  my $debT = $start+1; #sinon samtools faidx commence une base plus tôt (décalage)
	  print $fh "$chr\t$debT\t$endNe\t",distance($R1_uc,$R2_uc);
	  if ($motifPrint){
	    print $fh "\t$R1\t$R2";
	  }
	  print $fh "\n";
	}
      }
    }
  }
}

#aln("TTAACTTTT","aaaagtgaa");
sub aln {
  my ($s1,$s2)=@_;
#  print "S1 : $s1, S2 : $s2\n";
  my @parms;
  my $cpt=0;
  $s1=revCompMotif($s1);
  my ($align1, $align2);
  for (my $i = 1; $i <= length($s1); $i++){
    my $letter1 = substr($s1, $i-1, 1);
    my $letter2 = substr($s2, $i-1, 1);
    $letter2=uc($letter2);
    if ($letter1 ne $letter2){
      $cpt++;
      #print "$letter1\t$letter2\n";
    }
    last if $cpt==2;
    $align1 .= $letter1;
    $align2 .= $letter2;
  }
  #last if $cpt==2;
  #print "Alignement : \n$align1\n$align2";
  return ($align1,$align2);
}
