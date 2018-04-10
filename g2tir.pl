#!/usr/bin/perl
use strict;
use warnings 'all';
use Bio::SeqIO;
use Data::Dumper;
use Bio::DB::Fasta;
use Bio::Perl;
#use Devel::NYTProf;

use Getopt::Long qw/:config bundling auto_abbrev permute/;

# lancement/usage :
# perl test.pl motif genome
# ex : perl test.pl TTAAC[CAT][CT].[TCA] /home/lhelou/data/human_data/grch38-p10.fa

my $fa = "chr1_2_test.fa";
#my $fa = "/home/lhelou/Data/chr22.fa";
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



sub iupac {
  my ($self) = @_;
  print "self : $self\n";
  my $re = join '', map $IUPAC{ $_ }, split '', $self;
  #print $re,"\n";
  return $re;
}

# print $motifDirect,"\n";
my $momo = "TTAACHYNH";
my $motifDirect = iupac($momo);
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


    my @startTab;

      while ($seq =~ m/$listP/gio){
        #  my $motif = qr/$listP/;
        #  print "$listP\n$motif\n";
        $chr =~ s/(chr* ).*/$1/;
        #print $chr;
         my $pos=pos($seq); #pos en partant de la fin de notre séquence


         my $start = $pos - $motifSize; #$lengthDiv;#- 8;
         my $end = $pos;
         my $deb = $end+$min_length;
#         print "$chr\t$pos\t$start\t$pos\n";
         #my $loo = $end+$max_length;
         #my $loo = $max_length;
         #print $loo,"\n";
#push @chr, $chr unless $data{$chr};
#unless ($start ~~ @startTab){
        #push (@startTab, $start) unless $start ~~@startTab;#
        unless ($start ~~ @startTab){
          push (@startTab, $start);
          #print $start,"\n";
#           my $reverse = revCompMotif($listP);
# print $reverse,"\n";

          my $truc = substr($seq,$start,$motifSize);
          #print "startM : $truc\n";
          my $reverse = revCompMotif($truc);
          #print "$reverse\n";

          my $substr = substr($seq, $deb,$max_length);

          #my $degRev = motifSingle($reverse);
          #my $listRev = join ('|', @$degRev);
          #print $listN,"\n";

          #while ($substr =~m/$listN/gio){
          while ($substr  =~m/$listN/gio){

            #print "$substr\n$listRev\n\n";
            my $posNe = pos($substr);
            my $start2=$posNe-$motifSize+$deb;
            my $endNe = $posNe+$deb;

            my $revsub=substr($seq,$start2, $motifSize);


            my $PBLE = substr($seq,$start,($endNe-$start));


            my ($a,$b) = aln($truc,$revsub);
            #print "$truc\n$revsub\n\n";
            #
            # print "taille truc :", length($truc),"\n";
            # print "taille a :", length($a),"\n";
            if (length($truc)== length($a)){
              my $debT = $start+1; #sinon samtools faidx commence une base plus tôt (décalage)
                  #  print "R1 : $motif1\nR2 : $motif2Temp\nREV: $motif2\n\n";
                  #  print $a,"\n",$b,"\n\n";
            print "$chr\t$debT\t$endNe";
            if ($motifPrint){
              print "\t$truc\t$revsub";
            }
            print "\n";
            }
         }
        }
      }


    #  $Data::Dumper::Sortkeys = 1;
    #  print Dumper \%h;
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

