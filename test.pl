#!/usr/bin/perl
use strict;
use warnings 'all';
use Bio::SeqIO;
use Data::Dumper;
use Bio::DB::Fasta;
use Bio::Perl;
use Devel::NYTProf;

# lancement/usage :
# perl test.pl motif genome
# ex : perl test.pl TTAAC[CAT][CT].[TCA] /home/lhelou/data/human_data/grch38-p10.fa

################################EXTRACTMOTIF###########################
#if ($#ARGV != 1){
#  die "\n\tUsage : perl script FileName.fasta outputName\n\n";
#}

#open (my $fh, '>', $ARGV[1]); #ouverture/ecriture du rapport

############ RECHERCHE TIR 1
############ my $motifDirect = "TTAAC[CAT][CT].[TCA]"; #mon motif direct
############ my $motifIndComp = "[TGA].[AG][ATG]GTTAA";
############ FIN

############ RECHERCHE TIR 2
############my $motifDirect = "TTAAC[CAT][CGT]...."; #mon motif direct
############my $motifIndComp = "....[GCA][ATG]GTTAA";

############ RECHERCHE TIR 3
#my $motifDirect = "TTAAC[CAT][CAT][GTA]...."; #mon motif direct
#my $motifIndComp = "....[CAT][GTA][ATG]GTTAA";

############ RECHERCHE TIR 4
#my $motifDirect = "TT[GT]A."; #mon motif direct
#my $fa = "chr1_2_test.fa";
my $fa = "chr1_2_test.fa";

#  my $max_length = 30;
#  my $min_length = 2;
  #
my $max_length = 4000;
my $min_length = 43;
##VRAI
my $motifDirect = "TTAAC[CAT][CT].[TCA]";
#my $fa = '/home/lhelou/data/human_data/grch38-p10.fa';

#my $motifDirect = $ARGV[0];
#my $fa = $ARGV[1];

sub tailleMotif {
  #print "taille motif\n";
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


sub revCompMotif { # fonction permettant de générer le reverse comp du motif étudié
#print "revcom\n";
  my ($motif) = @_;# possible de reverse .{2} mais pas possible de reverse [AC]{2} => need amélioration
  my $rev = reverse($motifDirect);
  $rev =~ tr/ACGT[]/TGCA][/;
  $rev =~ s/\}(.*)\{\./\.{$1\}/;
  return $rev;
}

#print $motifDirect,"\n", $motifIndComp;

my $lengthDiv = tailleMotif($motifDirect)-1;

# sub findMotif {
#   my ($motif, $genome) = @_;
#   my $seqio_obj = Bio::SeqIO->new(-file => "$fa", #on recupere la sequence depuis notre fichier (arg 0)
#            -format => "fasta" );
#
#   my %h;
#   my $numID=1;
#
#   while ( my $seq_obj = $seqio_obj->next_seq ) {
#
#     my $chr = $seq_obj->display_id;
#     my $seq = $seq_obj->seq;
#     my $sensMotif;
#     my $lengthDiv = tailleMotif($motif)-1; #taille du motif - 1, utile pour calculer le start (8 pour le motif 1)
#
#     my $motifIndComp = revCompMotif($motif);
#     my @motifs = ($motif, $motifIndComp);
#
#     my ($start, $end, $subSeq);
#
#     foreach my $m (@motifs){
#
#       while ($seq =~ m/$m/gi){
#         $chr =~ s/(chr* ).*/$1/;
#         #récupérer position de début de notre motif
#         my $pos=pos($seq); #pos en partant de la fin de notre séquence
#         #print "position : $pos\n";
#         my $start = $pos - $lengthDiv;#- 8;
#         #pour motif long : $pos - 21
#         my $end = $pos;
#         my $subSeq = $seq_obj->subseq($start,$end); #récupération de la sous séquence représentant le motif
#         if ($m eq $motifIndComp){
#           $sensMotif = "-";
#         #DEBUG
#         #print "\n=====DEBUT=====\n",$subSeq,"\n",$m,"\n",$motifIndComp,"\n======\n";
#         #DEBUG
#         }
#         elsif ($m eq $motif){
#           $sensMotif ="+";
#         #DEBUG
#         #print "\n=====DEBUT SENS PLUS=====\n",$subSeq,"\n",$m,"\n",$motifDirect,"\n======\n";
#         #DEBUG
#         }
#         print "$numID\t$chr\t$start\t$end\t$sensMotif\t$subSeq\n";
#         $h{$numID}{"chr"}=$chr;
#         $h{$numID}{"start"}=$start;
#         $h{$numID}{"end"}=$end;
#         $h{$numID}{"strand"}=$sensMotif;
#         $h{$numID}{"subSeq"}=$subSeq;
#         $numID++;
#       }
#     }
#   }
#
#
#   #  print Dumper(\%h);
#   return \%h;
# }

my %fff = findMotif($motifDirect,$fa);
#print Dumper $fff;

# sub motifDegenere {
#   my $cpt = 0;
# #  my $lengthDiv = tailleMotif($motif)-1;
#   my ($pattern,$name)=@_;
#   my %patterns;# = ($pattern);
#
#   # my $patTemp = $pattern;
#   # my $cptMOT = 0;
#   # #$cptMOT{$1}++ while ($patTemp=~ /\[\w*\]/g);
#   #  while ($pattern=~m/\[\w*\]/g){
#   #    $cptMOT++;
#   #  }#$seq =~ m/$motif/gi
#   #
#   #
#   #  for (my $j=0; $j<$cptMOT; $j++){
#   #    my $cro = $pattern;
#   #    while ($cro=~m/\[\w*\]/g){
#   #      $pattern =~ s/\[\w*\]/\[\.\]/g;
#   #
#   #
#   #    }
#   #    print $pattern,"\n";
#   #  }
#   #
#   #  print $patTemp;
#   # print $cptMOT;
#   #$patTemp =~ s/\[\w*\]/\[\.\]/g;
#   #print $pattern,"\n",$patTemp,"\n";
#
#
#
# #  print length($pattern);
#   for (my $i=0; $i<length($pattern); $i++){
#     my $id = $cpt.$name;
#     $cpt++;
#     my $copy = $pattern;
#     my $ss = substr ($pattern,$i,1);
#
#
#
#     if ( ($ss ne "]" ) && ($ss ne "[")){
#       substr ($copy,$i,1) = ".";
#       $patterns{$id}{'motif'}=$copy;
# #      push @patterns, $copy;
#     ###$patterns{'motif'}=$copy;
#     }
#   }
#
#
#   print Dumper (\%patterns);
#   return %patterns;
#   #print @patterns;
#   #return @patterns;
# }

#Donc en gros, tu regardes si tu ne trouves pas de crocher ouvrant et fermant.
#Si le caractère que tu rencontres n'est ni ] ni [, alors du coup on remplace le caractère par un point.
#par contre si tu rencontres un crochet, tu incrémentes le compteur et là je comprends pas,
#tu veux faire quoi ? car là tu réafectes à la variable $ss une seule lettre.
#puis ensuite tu push le numéro de compteur au tableau pour ...?
# En gros je pige pas car vu le premier if, même après un "]" on va push nan ?


sub motifDegenere {
  #print "motifDegenere\n";
  my $cpt = 0;
  my ($pattern,$name)=@_;
  my %patterns;

  for (my $i=0; $i<length($pattern); $i++){
    my $id = $cpt.$name;
    $cpt++;
    my $copy = $pattern;
    my $ss = substr ($pattern,$i,1);

    if ( ($ss ne "]" ) && ($ss ne "[")){
      substr ($copy,$i,1) = ".";
      $patterns{$id}{'motif'}=$copy;
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
    $patterns{$id}{'motif'}=$copy;
    }
  }
#  print Dumper (\%patterns);
  return %patterns;
}


sub findMotif {
#  print "findMotif\n";
  my ($motif, $genome) = @_;
  my %sequences;

  my %hashChr;


  my $seqio_obj = Bio::SeqIO->new(-file => "$fa", #on recupere la sequence depuis notre fichier (arg 0)
           -format => "fasta" );

  my $motifIndComp = revCompMotif($motif);
#  my @motifs = ($motif, $motifIndComp);

  my $numID=1;
  my %motifPlusD = motifDegenere($motif,"p");
  foreach my $id (sort keys %motifPlusD){
    $motifPlusD{$id}{'sens'}="+";
  }
  my %motifMInusD = motifDegenere($motifIndComp,"n");
  foreach my $id (sort keys %motifMInusD){
    $motifMInusD{$id}{'sens'}="-";
  }
  my %allMotif = (%motifPlusD, %motifMInusD);


  while ( my $seq_obj = $seqio_obj->next_seq ) {

    my $chr = $seq_obj->display_id;
  #  print $chr,"\n";
    my $seq = $seq_obj->seq;

    my ($start, $end, $subSeq);



    foreach my $keysM (keys %allMotif){
    #  print "$keysM\n";

     my $motif = $allMotif{$keysM}->{motif};
     my $strand = $allMotif{$keysM}->{sens};


     #print "$motif\n";

     while ($seq =~ m/$motif/gi){
       #print "$motif\n";
       $chr =~ s/(chr* ).*/$1/;
       #récupérer position de début de notre motif
       my $pos=pos($seq); #pos en partant de la fin de notre séquence
       #print "position : $pos\n";
       my $start = $pos - $lengthDiv;#- 8;
       #pour motif long : $pos - 21
       my $end = $pos;
       ##### my $subSeq = $seq_obj->subseq($start,$end); #récupération de la sous séquence représentant le motif
       ###DEBUG print "$subSeq\n\n";

       ##CONSTRUCTION table de hach contenant tous les motifs
       #my $id = $cpt;
       ##print $id,"\n";

          $hashChr{$chr}{$start}{'end'}=$end;
          $hashChr{$chr}{$start}{'strand'}=$strand;

      }
    }
  }
  #print Dumper \%allMotif;
 print Dumper \%hashChr;
  return %hashChr;

}

##open (my $fh, '>', "test.out");

#Fonction qui paire les motifs
#Première chose à faire est de récupérer le start que si le sens est + et le end si le sens est -
sub pairsMotif {
#  print "pairs motif\n";
  my (%h)=@_;
  my (%data, @chr, @tab);

  #DEBUG
#  print Dumper \%h;

   foreach my $chrH (sort keys %h){
  # #print $id,"\n";
   #print $startK,"\n";
   my $chr = $chrH;
   foreach my $stH (sort keys %{$h{$chrH}}){
     #print $stH,"\n";

   my $start = $stH;
   my $end = $h{$chrH}{$stH}->{end};
   my $strand = $h{$chrH}{$stH}->{strand};
#   print $start,"\t",$end,"\n";
  # my $end = ${}
   push @chr, $chr unless $data{$chr};
   push @{$data{$chr}}, $strand eq '+' ? [$start, $strand] : [$end, $strand];
   }
 }

  #foreach my $t (keys %{%$h{$id}}){
  #  print "ttttttttttttt", $t,"\n";
  # my $chr = $startK;
  # my $strand = %$h{$startK}->{strand};
  # my $start =  $startK;
  # my $end =  %$h{$startK}->{end};
  # my $subSeq =  %$h{$startK}->{subSeq};




    # my $chr = %$h{$startK}->{chr};
    # my $strand = %$h{$startK}->{strand};
    # my $start =  $startK;
    # my $end =  %$h{$startK}->{end};
    # my $subSeq =  %$h{$startK}->{subSeq};
    ####DEBUG
    #print %$h{$id}->{chr},"\t";
    #print %$h{$id}->{start},"\t";
    #print %$h{$id}->{strand},"\n\n";


    #
    # push @chr, $chr unless $data{$chr};
    # push @{$data{$chr}}, $strand eq '+' ? [$start, $strand] : [$end, $strand];
  #}

  for my $chr(@chr){
    my $data = $data{$chr};
    @$data = sort { $a->[0] <=> $b->[0] } @$data; #car séparé selon + et moins

    for my $i ( 0 .. $#$data ) { # $#$data est un élément de moins que dans @$data
    #print "le i : $i\n";

    next unless $data->[$i][1] eq '+';

    for my $j ( $i + 1 .. $#$data ) {
      next unless $data->[$j][1] eq '-';
      my $length = $data->[$j][0] - $data->[$i][0]; #$j end (strand -) et $i start strand +
      last if $length > $max_length;
        if ($length>$min_length){
        ### pour push de tableau dans un tableau cf programmation perl 3ème édition p249
          push @tab, [$chr, $data->[$i][0], $data->[$j][0]];
        }
      }
    }
  }
  # parcourir chaque/case ligne du tableau
  # for my $ligne (@tab){
  #   ##print $fh "@$ligne\n";
  #   print "@$ligne\n";
  # }
  return @tab;
}

############### BONNE FONCTION DEBUT
# sub pairsMotif {
# print "pairs motif\n";
#   my ($h)=@_;
#   my (%data, @chr, @tab);
#
#   #DEBUG
#   #print Dumper ($h);
#
#   foreach my $startK (sort keys %$h){
#   #print $id,"\n";
# print $startK,"\n";
#   foreach my $st (keys %{%$h{$startK}}){
#     print $st;
#   }
#
#   #foreach my $t (keys %{%$h{$id}}){
#   #  print "ttttttttttttt", $t,"\n";
#   my $chr = $startK;
#   my $strand = %$h{$startK}->{strand};
#   my $start =  $startK;
#   my $end =  %$h{$startK}->{end};
#   my $subSeq =  %$h{$startK}->{subSeq};
#     # my $chr = %$h{$startK}->{chr};
#     # my $strand = %$h{$startK}->{strand};
#     # my $start =  $startK;
#     # my $end =  %$h{$startK}->{end};
#     # my $subSeq =  %$h{$startK}->{subSeq};
#     ####DEBUG
#     #print %$h{$id}->{chr},"\t";
#     #print %$h{$id}->{start},"\t";
#     #print %$h{$id}->{strand},"\n\n";
#     push @chr, $chr unless $data{$chr};
#     push @{$data{$chr}}, $strand eq '+' ? [$start, $strand] : [$end, $strand];
#   }
#
#   for my $chr(@chr){
#     my $data = $data{$chr};
#     @$data = sort { $a->[0] <=> $b->[0] } @$data; #car séparé selon + et moins
#
#     for my $i ( 0 .. $#$data ) { # $#$data est un élément de moins que dans @$data
#     #print "le i : $i\n";
#
#     next unless $data->[$i][1] eq '+';
#
#     for my $j ( $i + 1 .. $#$data ) {
#       next unless $data->[$j][1] eq '-';
#       my $length = $data->[$j][0] - $data->[$i][0]; #$j end (strand -) et $i start strand +
#       last if $length > $max_length;
#         if ($length>$min_length){
#         ### pour push de tableau dans un tableau cf programmation perl 3ème édition p249
#           push @tab, [$chr, $data->[$i][0], $data->[$j][0]];
#         }
#       }
#     }
#   }
#   # parcourir chaque/case ligne du tableau
#   # for my $ligne (@tab){
#   #   ##print $fh "@$ligne\n";
#   #   print "@$ligne\n";
#   # }
#   return @tab;
# }

###################### BONNE FONCTION FIN




# sub pairsMotif {
#
#   my ($h)=@_;
#   my (%data, @chr, @tab);
#
#   #DEBUG
#   #print Dumper ($h);
#
#   foreach my $id (sort keys %$h){
#   #print $id,"\n";
#
#   #foreach my $t (keys %{%$h{$id}}){
#   #  print "ttttttttttttt", $t,"\n";
#     my $chr = %$h{$id}->{chr};
#     my $strand = %$h{$id}->{strand};
#     my $start =  %$h{$id}->{start};
#     my $end =  %$h{$id}->{end};
#     my $subSeq =  %$h{$id}->{subSeq};
#     ####DEBUG
#     #print %$h{$id}->{chr},"\t";
#     #print %$h{$id}->{start},"\t";
#     #print %$h{$id}->{strand},"\n\n";
#     push @chr, $chr unless $data{$chr};
#     push @{$data{$chr}}, $strand eq '+' ? [$start, $strand] : [$end, $strand];
#   }
#
#   for my $chr(@chr){
#     my $data = $data{$chr};
#     @$data = sort { $a->[0] <=> $b->[0] } @$data; #car séparé selon + et moins
#
#     for my $i ( 0 .. $#$data ) { # $#$data est un élément de moins que dans @$data
#     #print "le i : $i\n";
#
#     next unless $data->[$i][1] eq '+';
#
#     for my $j ( $i + 1 .. $#$data ) {
#       next unless $data->[$j][1] eq '-';
#       my $length = $data->[$j][0] - $data->[$i][0]; #$j end (strand -) et $i start strand +
#       last if $length > $max_length;
#         if ($length>$min_length){
#         ### pour push de tableau dans un tableau cf programmation perl 3ème édition p249
#           push @tab, [$chr, $data->[$i][0], $data->[$j][0]];
#         }
#       }
#     }
#   }
#   # parcourir chaque/case ligne du tableau
#   for my $ligne (@tab){
#     ##print $fh "@$ligne\n";
#     print "@$ligne\n";
#   }
#   return @tab;
# }

my @tab = pairsMotif(%fff);
#__END__
#filterMotif(\@tab);

sub filterMotif {
  print "filtre motif\n";
  #print $lengthDiv,"\n";
  my ($in)=@_;
  my $db = Bio::DB::Fasta -> new($fa);
  my @tabou;
  for my $ligne (@$in){
  #print "@$ligne\n";
    my $chr = @$ligne[0];
    my $start = @$ligne[1];
    my $end = @$ligne[2];
    #    my ($chr, $start,$end)=split " ", @$ligne;
    #print "$chr\t$start\t$end\n";
    my ($end1,$start2)=0;
    $end1 = $start+$lengthDiv;
    $start2 = $end-$lengthDiv;
    my $motif1 = $db->subseq($chr, $start,$end1);
    my $motif2Temp = $db->subseq($chr, $start2,$end);
    my $motif2 = revcom($motif2Temp) -> seq;

    my ($a,$b) = aln($motif1,$motif2);

    if (length($motif1)== length($a)){
          #  print "R1 : $motif1\nR2 : $motif2Temp\nREV: $motif2\n\n";
          #  print $a,"\n",$b,"\n\n";
           push @tabou, [$chr, $start, $end, $motif1, $motif2Temp];
    }
  }

    for my $ligne (@tabou){
    ##print $fh "@$ligne\n";
    print "@$ligne\n";
  }
    return @tabou;
}

sub aln {
  my ($s1,$s2)=@_;
  my @parms;
  my $cpt=0;
  my ($align1, $align2);
  for (my $i = 1; $i <= length($s1); $i++){
    my $letter1 = substr($s1, $i-1, 1);
    my $letter2 = substr($s2, $i-1, 1);
    if ($letter1 ne $letter2){
      $cpt++;
      #print "$letter1\t$letter2\n";
    }
    last if $cpt==2;
    $align1 .= $letter1;
    $align2 .= $letter2;
  }
  #last if $cpt==2;
  return ($align1,$align2);
}
