#!/usr/bin/perl
use strict;
use warnings 'all';
use Bio::SeqIO;
use Data::Dumper;
use Bio::DB::Fasta;
use Bio::Perl;
use Text::Levenshtein::XS qw/distance/;
use Getopt::Long qw/:config bundling auto_abbrev permute/;
use Pod::Usage;

# lancement/usage :
# perl g2tir.pl --pm --motif motif_iupac --in fasta_file

# time perl ~/scripts/perl/g2TIR/g2tir/g2tir.pl --score 0 --min 800  --printScore --config --pm --motif TTAANNNN --in /home/lhelou/ancestralSpecies/Analyses/phylogenieSequencesProches/pgbd5_vs_uro/outDir/debug_g2tir/QLIY01000744_1_57307_64203.fa --out t5

my $max_length = 6000;
my $min_length = 43;

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


my ($motifPrint,$motifDegenerated, $iupac, $fa, $out_fasta, $sum, $fhfasta, $printScore, $TIR1, $TIR2, $score, $TSD);

my $time = time;

my $output = "g2tir_out.$time";

my ($nomprog) = $0 =~ /([^\/]+)\d*$/;
my $cmd_line = $nomprog." @ARGV";

my $help = 0;

GetOptions ("help|?"=> \$help,
            "in=s"=> \$fa,
            "motif=s"=> \$iupac,
            "out=s"=> \$output,
            "maxL=i"=> \$max_length,
            "minL=i"=> \$min_length,
            "pm" => \$motifPrint, #print les TIR
            "degen" => \$motifDegenerated, #dégénère le motif
            "fasta|fa"=> \$out_fasta,
            "config"=> \$sum, #file output config
            "printScore"=> \$printScore,
            "TIR1|R1=s"=> \$TIR1,
            "TIR2|R2=s"=> \$TIR2,
            "score=i"=> \$score,
            "TSD=i"=> \$TSD #print les bases à l'extérieur
            )
            or pod2usage(2);
pod2usage(1) if $help;

pod2usage(2) if ! defined ($fa && $iupac || $TIR1);

# die "You must specify a fasta file (option --in fasta_file)\n" if ! defined $fa;
# die "You must specify a motif in iupac format (option --motif iupac_motif)\n" if ! defined ($iupac || $TIR1);
# die "The fasta file <$fa> doesn't exist" if ! -e $fa;
# die "Maximum and minimum lengths are wrong, maxL must be higher than minL\n" if $max_length<$min_length;

# die "You must specify a fasta file (option --in fasta_file)\n" if ! defined $fa;
# die "You must specify a motif in iupac format (option --motif iupac_motif)\n" if ! defined ($iupac || $TIR1);
# die "The fasta file <$fa> doesn't exist" if ! -e $fa;
# die "Maximum and minimum lengths are wrong, maxL must be higher than minL\n" if $max_length<$min_length;

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

# if ($TIR1){
#   $motifDirect = iupacToRegex($R1);
# }
# elsif ($TIR1){
#   $motifDirect = iupacToRegex($TIR1);
#   $motifIndirect = iupacToRegex($TIR2);
# }

if ($iupac){
  $motifDirect = iupacToRegex($iupac);
}
elsif ($TIR1){
  $motifDirect = iupacToRegex($TIR1);
  $motifIndirect = iupacToRegex($TIR2);
}

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
    $revmo = $motifIndirect;
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



my $try = findPairsMotif($tabP, $tabN);

sub findPairsMotif {
  my ($tabP, $tabN)= @_;
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
	           last;
          }
        }
        @endPos = $index < scalar @endPos ? @endPos[$index..$#endPos] : ();
      }
      $minSkyline = $pEnd + $min_length;
      my $R1 = substr($seq, $pStart, $motifSizeDirect);
  #    print "R1 :",$R1,"\n";
      my $TSD_R1;
      if ($TSD) {
        $TSD_R1 = substr($seq,$pStart-$TSD,$TSD);
      }
      my $R1_uc = uc($R1);
      if (@endPos){
        if ($startSubstr < length($seq)){
      	$startSubstr = $endPos[-1] ;
	}
      	$length_extension = $max_length - ($startSubstr - $pEnd - $min_length);
      }
      else {
      	$startSubstr = $minSkyline;

      	$length_extension = $max_length - $min_length;
}
	my $substr;
	next if ($startSubstr > length($seq));

	if ($startSubstr + $length_extension > length($seq)){
		$length_extension = length($seq)-$startSubstr;
	}

      $substr = substr($seq, $startSubstr, $length_extension);

      while ($substr =~ m/(?=($listN))/gio){#} && defined $substr){
  ##TRUE    while ($substr =~ m/$listN/gio){#} && defined $substr){
        my $nPos=pos($substr)+$motifSizeDirect;
##TRUE        my $nPos = pos($substr);
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

        my $TSD_R2;
        if ($TSD) {
#          $TSD_R2 = substr($seq,$nStart+$motifSizeIndirect+$TSD,$TSD);
          $TSD_R2 = substr($seq,$nEndPos,$TSD);
        }

      	my $R2_uc = revCompMotif(uc($R2));

	# if ($R1_uc eq $R2_uc){
	# 	print "R1 : $R1_uc\nR2 : $R2_uc\nYOUHOUUU\n";
	# }
      	my $dist = distance($R1_uc, $R2_uc);

        my $scoreComp;#=0;

        if (defined $score){
          $scoreComp = $score;
        }
        else {
          $scoreComp = $maxScore;
          if ($R1){
            $scoreComp = 200;
          }
        }

    #    if ($iupac){
            if ($dist <= $scoreComp){
    #    my @out = ($chr, $pStart, $nEndPos, $nEndPos- $pStart, $dist);
    #          $pStart=$pStart+1; #mise à niveau pour standart extraction fasta
              my @out = ($chr, $pStart, $nEndPos);
              push @out, ($R1, $R2) if $motifPrint;
              push @out, ($dist) if $printScore;
              push @out, ($TSD_R1, $TSD_R2) if $TSD;
              print $fh join ("\t", @out), "\n";
              if ($out_fasta){
    #            $pStart =$pStart+1;
                my $PBLE = $seq_obj->subseq($pStart+1, $nEndPos);
                print $fhfasta ">$chr:$pStart-$nEndPos\n$PBLE\n";
              }
            }
          #}

  #       elsif ($TIR1){
  #         if ($dist <= $scoreComp){
  #           my @out = ($chr, $pStart, $nEndPos);
  #           push @out, ($R1, $R2) if $motifPrint;
  #           push @out, ($dist) if $printScore;
  #           push @out, ($TSD_R1, $TSD_R2) if $TSD;
  #           print $fh join ("\t", @out), "\n";
  #           if ($out_fasta){
  # #            $pStart =$pStart+1;
  #             my $PBLE = $seq_obj->subseq($pStart+1, $nEndPos);
  #             print $fhfasta ">$chr:$pStart-$nEndPos\n$PBLE\n";
  #           }
  #         }
  #       }
#}

      }
    }
  }
}

if ($sum){
  my $file_sum = "$output.config";
  open (my $fhsum, '>', $file_sum);
  print $fhsum "cmd line : $cmd_line\n";
  print $fhsum "motif : $iupac\nfasta : $fa\n" if ($iupac);
  print $fhsum "motif was degenerated\n" if ($motifDegenerated);
  print $fhsum "fasta file was generated\n" if ($out_fasta);
}

__END__

=head1 G2TIR

sample - Using truc truc

=head1 SYNOPSIS

g2tir.pl --score 0 --motif --in file.fasta --out output_file

--help for more options

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exits.

=item B<--in> Fasta file [required]

Fasta file in which the ITRs are searched for

=item B<--motif> motif for ITR in iupac format [required]

=item B<--R1> motif for ITR1 in iupac format if R1 and R2 are different [required]

=item B<--R2> motif for ITR2 in iupac format if R1 and R2 are different [required]

=item B<--score>

=item B<--out> output file with the results in table format

=item B<--maxL> Maximum distance between two ITRs (INT default : 6400 nt)

=item B<--minL> Minimum distance between two ITRs (INT default : 43 nt)

=item B<--pm> Print motifs of ITR founded in the output file

=item B<--printScore> Print the score in the output file

=item B<--TSD> Print the TSD (4nt before each TIR) in the output file

=item B<--degen> Degenerate the motif submitted

=item B<--fasta> Generate a fasta file with the sequences identified

=item B<--config> Generate a config file with the option selected

=cut

#time perl ~/scripts/perl/g2TIR/g2tir/g2tir_impr.pl --printScore --fasta --config --pm --motif TTAANNNN --in sra_tblastn_wgs.fasta --score 0 --out sra_tblastn_wgs/g2tir_results/TTAANNNN &
