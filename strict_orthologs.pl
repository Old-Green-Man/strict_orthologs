#!/usr/bin/env perl

my @BOOL = qw( false true );
my %opts = ( DEBUG => 0, VERBOSE => 0, METHOD => 'complete_linkage', USE_PARALOGS => 1, MIN_BITS => 40, MIN_PCNT => 0, MIN_PCNT2 => 0, DIFF_PCNT => 20, MIN_ALEN => 10, SINGLETS => 1, SEP => "_", PARSED_HITS => undef, PARSED_SELF => undef, MIN_PARALOG_PCNT => 50.0 );

use warnings;
use strict;
use Getopt::Long;
use Time::HiRes qw( time );

my $OLD_WAY = 0;

srand(11);
my $t_start = time;
my ( $seqLenFile, @tabularBlastFiles ) = parseArgs();
my $o = DnaNinja::Ortho->new( \%opts );

my $start = time;
$o->getSeqLengths( $seqLenFile );
$o->logDurationMemory( "getSeqLengths", $o->sec2time(time - $start) );

if( defined $o->{PARSED_HITS} ) {
  $start = time;
  $o->parsePrecompiledData( $o->{PARSED_HITS}, $o->{PARSED_SELF} );
  $o->logDurationMemory( "parsePrecompiledData", $o->sec2time(time - $start) );
} elsif( defined $o->{RBH_FILE} ) {
  $start = time;
  $o->parseRecipricalBestHits( $o->{RBH_FILE} );
  $o->logDurationMemory( "parseRecipricalBestHits", $o->sec2time(time - $start) );
} else {
  $start = time;
  foreach my $file ( @tabularBlastFiles ) {
    $o->parseTabularBlast( $file );
  }
  $o->logDurationMemory( "parseTabularBlast", $o->sec2time(time - $start) );
}
$start = time;
if( $o->{USE_PARALOGS} ) {
  if( $OLD_WAY ) {
    $o->groupParalogs_old;
    $o->logDurationMemory( "groupParalogs_old", $o->sec2time(time - $start) );
  } else {
    if( not $o->{RBH_FILE} ) {
      $o->groupParalogs;
      $o->logDurationMemory( "groupParalogs", $o->sec2time(time - $start) );
      $start = time;
    }
    $o->findOrthos;
    $o->logDurationMemory( "findOrthos", $o->sec2time(time - $start) );
  }
} else {
  $o->findOrthos;
  $o->logDurationMemory( "findOrthos", $o->sec2time(time - $start) );
}

if( not defined $o->{RBH_FILE} ) {
  $start = time;
  $o->writeRecipricalBestHits;
  $o->logDurationMemory( "writeHitsFile", $o->sec2time(time - $start) );
}

$start = time;
$o->clusterCompleteLinkage();
$o->logDurationMemory( "clusterCompleteLinkage", $o->sec2time(time - $start) );

$start = time;
$o->printClusters();
$o->logDurationMemory( "printClusters", $o->sec2time(time - $start) );

$o->logDurationMemory( "Total", $o->sec2time(time - $t_start) );

sub parseArgs {
  my @args = qw( gene-lengths-file );
  my $usage = qq{
usage: $0 @args [tabular-blast-output-file tabular-blast-output-file2 ...] [options]

This reads tabular blast output. It assumes that blast was run with the
following -outfmt:
-outfmt="6 qacc sacc bitscore nident qstart qend sstart send"

When run the first time on blastoutput it will write a file rbh.dat.  This can be read in using -rbh_file=rbh.dat without needing the blast outputs and speeding things up.

gene-lengths-file should contain the lengths of all genes (format: len<TAB>id)

options
-------
-rbh_file preparsed reciprical best hits including paralog links.  This bypasses
 running groupParalogs and should include all relevant reciprical best hits
 needed for clustering.

-parsedHits preparsed %HITS data - deprecated

-parsedSelf preparsed %PARALOGS data - deprecated

-sep (default: $opts{SEP}) separator between taxon name and gene in the query
 and sbjct ids. For example, genomeA_geneX.

-useParalogs (default: $BOOL[$opts{USE_PARALOGS}] Include paralogs that are specific to a
 taxon into the same ortholog group.  These paralogs will be considered a single
 gene internally that share the various reciprical hits.

-minParalogPcnt (default: $opts{MIN_PARALOG_PCNT}) if n_ident / max(q,s) * 100 is less than this then not counted as paralogs

-singlets (default: $BOOL[$opts{SINGLETS}]) print singlets as their own group

-minPcnt  (default: $opts{MIN_PCNT}) minimum value of 100 * identities / max( qlen, slen ) to keep
  NOT CURRENTLY USED

-minPcnt2 (default: $opts{MIN_PCNT2}) minimum value of 100 * identities / qcov to keep
 This second (normal percent identity) is used linking CRISPR arrays where we
 don't want penalize so much for subsequence matches as long as a sufficient
 amount is aligned as determined by -minAlen
  NOT CURRENTLY USED

-minBits  (default: $opts{MIN_BITS}) minimum bit value to keep

-minAlen (default: $opts{MIN_ALEN}) minimum per HSP alignment length.  This is to
 capture alignments of adjacent repeats and spacers.  The length should be at least
 greater than the maximum expected length of either repeat or spacer.

-diffPcnt (default: $opts{DIFF_PCNT}) if an HSP differs in percent identity from the
 first HSP for a given hit by more than this exclude it from calculating overall
 pcntId, bits, and coverage.

-help

-verbose (default: $BOOL[$opts{VERBOSE}]) verbose logging

-debug   (default: $BOOL[$opts{DEBUG}]) Run extra debug checks

};

  # -method not implemented yet - 'single_linkage', 'complete_linkage', 'best_score' (default: $opts{METHOD})
  # best_score keeps a single gene from a genome in cluster based on alignment score

  my $help = 0;
  my $result = GetOptions
    (
     'rbh_file=s' => \$opts{RBH_FILE},
     'parsedHits=s' => \$opts{PARSED_HITS},
     'parsedSelf=s' => \$opts{PARSED_SELF},
     'sep=s' => \$opts{SEP},
     'singlets!' => \$opts{SINGLETS},
     'useParalogs!' => \$opts{USE_PARALOGS},
     'minParalogPcnt=f' => \$opts{MIN_PARALOG_PCNT},
     'minPcnt=f' => \$opts{MIN_PCNT},
     'minPcnt2=f' => \$opts{MIN_PCNT2},
     'minBits=f' => \$opts{MIN_BITS},
     'minAlen=i' => \$opts{MIN_ALEN},
     'diffPcnt=f' => \$opts{DIFF_PCNT},
     'help!' => \$opts{help},
     'verbose!' => \$opts{VERBOSE},
     'debug!' => \$opts{DEBUG},
     'method=s' => \$opts{METHOD},
    );
  $result or print $usage and exit;
  $help and print $usage and exit;
  if( defined $opts{PARSED_HITS} or defined $opts{PARSED_SELF} ) {
    if( not defined $opts{PARSED_SELF} and $opts{USE_PARALOGS} ) {
      die "error: if -parsedHits is used and using paralogs then -parseSelf must be specified\n";
    }
    @ARGV == 1 or die "error: if -parsedHits is used then there should be only 1 non option argument. However, we should allow mixing of -parsedHits and blast outputs\n";
  } elsif( defined $opts{RBH_FILE} ) {
    @ARGV == 1 or die "error: if -rbh_file is used then the only argument should be the sequence lengths file. However, we should later allow mixing of -rbh_file and blast outputs\n";
  } else {
    @ARGV >= @args or print $usage and exit;
  }
  $opts{MIN_PCNT} /= 100;
  $opts{MIN_PCNT2} /= 100;
  $opts{DIFF_PCNT} /= 100;
  return @ARGV;
}

package DnaNinja::Ortho;

use warnings;
use strict;
use Time::HiRes qw( time );
use Carp;

# Some options are required by certain methods. %opts can passed to new to set these options.  For example:
# %opts = ( DEBUG => 0, VERBOSE => 0, METHOD => 'complete_linkage', USE_PARALOGS => 1, MIN_BITS => 40, MIN_PCNT => 0, MIN_PCNT2 => 0, DIFF_PCNT => 20, MIN_ALEN => 10, SINGLETS => 1, SEP => "_", PARSED_HITS => undef, PARSED_SELF => undef, MIN_PARALOG_PCNT => 50.0 );
# I need to give these default values are require them to be passed the constructor.
sub new {
  my $class_or_object = shift;
  my $opts = shift;
  my $class = ref($class_or_object) || $class_or_object;
  my $self = {};
  bless $self, $class;
  $self->_initialize( $opts );
  return $self;
}

sub _initialize {
  my $o = shift;
  my $opts = shift;
  foreach my $opt ( keys %$opts ) {
    $o->{$opt} = $opts->{$opt};
  }
}

sub logDurationMemory {
  my $o = shift;
  my $str = shift;
  my $time = shift;
  printf STDERR "Duration: $time";
  my $mem = $o->getMemory;
  defined $mem and printf STDERR "\tMemory: $mem";
  print STDERR "\t$str\n";
  STDERR->flush;
}

sub getMemory {
  my $o = shift;
  open( STAT , "</proc/$$/status" ) or die "Unable to open stat file";
  my ( @buf );
  while( <STAT> ) {
    /^VmSize:\s*(\d+)/ and push @buf, sprintf "VmSize: %.2f", $1 / 1024**2;
    if( /^VmRSS:\s*(\d+)/ ) {
      push @buf, sprintf "VmRSS: %.2f", $1 / 1024**2;
      last;
    }
  }
  @buf or return undef;
  return sprintf "%s", join " ", @buf;
}

sub timeStamp {
  my $o = shift;
  my $str = shift;
  printf STDERR "%s $str\n", scalar localtime;
  STDERR->flush;
}

sub sec2time {
  my $o = shift;
  my $sec = shift;
  my ( $hours, $min );
  if( $sec / 3600 >= 1 ) {
    $hours = int $sec / 3600;
    $sec = $sec % 3600;
    $min = int $sec / 60;
    $sec = $sec % 60;
    return sprintf "%02d:%02d:%.4g", $hours, $min, $sec;
  } elsif ( $sec / 60 >= 1 ) {
    $min = int $sec / 60;
    $sec = $sec % 60;
    return sprintf "00:%02d:%.4g", $min, $sec;
  } else {
    return sprintf "00:00:%.4g", $sec;
  }
}

sub parseRecipricalBestHits {
  my $o = shift;
  my $rbh_file = shift;
  my $fh = getFH( $rbh_file );
  my ($qid, $sid, $bits, $pcnt);
  while( <$fh> ) {
    chomp;
    ($qid, $sid, $bits, $pcnt) = split /\t/;
    $o->{BEST_HIT}{$qid}{$sid} = [ $bits, $pcnt ];
  }
}

sub parsePrecompiledData {
  my $o = shift;
  my $hitsfile = shift;
  my $selffile = shift;
  my $fh = getFH( $hitsfile );
  my ($qid, $sid, $bits, $pcnt);
  while( <$fh> ) {
    chomp;
    ($qid, $sid, $bits, $pcnt) = split /\t/;
    $o->{BEST_HIT}{$qid}{$sid} = [ $bits, $pcnt ];
  }
  defined $selffile or return;
  $fh = getFH( $selffile );
  while( <$fh> ) {
    chomp;
    ($qid, $sid, $bits, $pcnt) = split /\t/;
    $o->{PARALOGS}{$qid}{$sid} = [ $bits, $pcnt ];
  }
}

# This assumes all data from a particular genome vs another genome are together,
# but doesn't assume anything about the order of hits or hsps within that
# group. We require the genome vs genome data to be grouped so that we can pipe
# in large numbers of blast searches but not have to save it all in memory.
# This assumes tabular blast fields were generated with the following -outfmt="6
# qseqid sseqid bitscore nident qstart qend sstart send"
sub parseTabularBlast {
  my $o = shift;
  my $file = shift;
  my $fh = getFH( $file );
  my ( %hsps );
  $_ = <$fh>; s/\s+$//; # first line
  my ( $qid, $sid, $bits, $nident, $qb, $qe, $sb, $se ) = split /\t/;
  my $qid_last = $qid;
  push @{$hsps{$sid}}, [ $bits, $nident, $qb, $qe, $sb, $se ];
  while ( <$fh> ) {
    s/\s+$//;
    ( $qid, $sid, $bits, $nident, $qb, $qe, $sb, $se ) = split /\t/;
    if( $qid ne $qid_last ) {
      $o->processQuery( $qid_last, \%hsps );
      $qid_last = $qid;
      %hsps = ();
    }
    push @{$hsps{$sid}}, [ $bits, $nident, $qb, $qe, $sb, $se ];
  }
  $o->processQuery( $qid_last, \%hsps ); # end condition
}

# This will usually be one genome_gene vs some genome, but could, in principle, be
# one genome_gene vs multiple genomes.
sub processQuery {
  my $o = shift;
  my $qid = shift;
  my $hsps = shift;
  my ( $q_taxon, $q_gene ) = $qid =~ /(\S+)$o->{SEP}(\S+)/;
  my ( $bits, $pcnt, $nident, $alen, $nhsps );
  my ( $s_taxon, $s_gene, $sid, %best_hit );
  foreach $sid ( sort keys %$hsps ) { # sort is just to make self.dat output reproducible
    ( $s_taxon, $s_gene ) = $sid =~ /(\S+)$o->{SEP}(\S+)/;
    $nhsps = @{$hsps->{$sid}};
    if( $nhsps == 1 ) {
      ( $bits, $nident ) = ( $hsps->{$sid}[0][0], $hsps->{$sid}[0][1] );
    } else {
      ( $bits, $nident, undef, undef ) = mergeHsps( $hsps->{$sid}, { diff_pcnt => $o->{DIFF_PCNT} } );
    }
    defined $bits or next;
    $bits < $o->{MIN_BITS} and next;
    # The max( $o->{LEN}{$qid}, $o->{LEN}{$sid} ) here is to penalize proteins/paralogs that differ greatly in length.
    $pcnt = sprintf "%.4g", 100 * $nident / max( $o->{LEN}{$qid}, $o->{LEN}{$sid} );
    if( $q_taxon eq $s_taxon ) {
      $q_gene eq $s_gene and next; # self hit
      $pcnt < $o->{MIN_PARALOG_PCNT} and next; # don't count distant homologs as close paralogs
      $o->{PARALOGS}{$qid}{$sid} = [ $bits, $pcnt ];
      # print {$o->{fh_self}} "$qid\t$sid\t$bits\t$pcnt\n";
    } else {
      # This keeps only 1 best hit per taxon. Make sure this doesn't mess up paralog stuff.
      if( defined $best_hit{$s_taxon} ) {
        if($bits <= $best_hit{$s_taxon}[1]) {
          # not best hit
          # If it's equal it's a close paralog to best and will be caught when doing paralogs
          next;
        }
        $o->{BEST_HIT}{$qid}{$sid} = [ $bits, $pcnt ];
        delete $o->{BEST_HIT}{$qid}{$best_hit{$s_taxon}[0]};
        $best_hit{$s_taxon} = [$sid, $bits];
      } else {
        $o->{BEST_HIT}{$qid}{$sid} = [ $bits, $pcnt ];
        $best_hit{$s_taxon} = [$sid, $bits];
      }
    }
  }
}

sub clusterCompleteLinkage {
  my $o = shift;
  my ($ortho, $idA, $idB, $bits, $removed, @removed, $start );
  my $cnt = 0;

  while ( $ortho = shift @{$o->{ORTHOS}} ) {
    ($bits, $idA, $idB) = @$ortho;

    if (exists $o->{ID2CLSTR}{$idA} and exists $o->{ID2CLSTR}{$idB}) { # both genes in SOME cluster due to links to other proteins
      if ($o->{ID2CLSTR}{$idA} != $o->{ID2CLSTR}{$idB}) {
        # In different clusters.  See if they can be merged.
        exists $o->{nomerge}{$o->{ID2CLSTR}{$idA}}{$o->{ID2CLSTR}{$idB}} and next;
        $start = time;
        $removed = $o->checkCompleteMerge($o->{ID2CLSTR}{$idA}, $o->{ID2CLSTR}{$idB});
        $o->{time}{checkCompleteMerge} += time - $start;
        $o->{calls}{checkCompleteMerge}++;
        defined $removed and push @removed, $removed;
      }
    } elsif (exists $o->{ID2CLSTR}{$idA}) { # gene $idB not in cluster
      $start = time;
      $o->checkCompleteLinkage( $idB, $o->{ID2CLSTR}{$idA} );
      $o->{time}{checkCompleteLinkage} += time - $start;
        $o->{calls}{checkCompleteLinkage}++;
    } elsif (exists $o->{ID2CLSTR}{$idB}) { # gene $idA not in cluster
      $start = time;
      $o->checkCompleteLinkage( $idA, $o->{ID2CLSTR}{$idB} );
      $o->{time}{checkCompleteLinkage} += time - $start;
      $o->{calls}{checkCompleteLinkage}++;
    } else { # neither of the genes, which are reciprical best hits, are in cluster, create new cluster
      $start = time;
      $o->{VERBOSE} and printf STDERR "INFO: clusterCompleteLinkage creating oid=%d with $idA and $idB\n", $cnt;
      $o->{ID2CLSTR}{$idA} = $cnt;
      $o->{ID2CLSTR}{$idB} = $cnt;
      $o->{CLUSTERS}[$cnt]{links}{$idA}{$idB} = undef;
      $o->{CLUSTERS}[$cnt]{links}{$idB}{$idA} = undef;
      $cnt++;
      # while( defined $o->{CLUSTERS}[$cnt] ) { $cnt++ }
      $o->{time}{new_cluster} += time - $start;
      $o->{calls}{new_cluster}++;
    }
  }
  foreach $cnt ( sort { $b <=> $a } @removed ) {
    defined $o->{CLUSTERS}[$cnt] and die "error: \$o->{CLUSTERS}[$cnt] shouldn't be defined\n";
    splice @{$o->{CLUSTERS}}, $cnt, 1;
  }
  foreach my $func ( qw(checkCompleteMerge checkCompleteLinkage new_cluster) ) {
    printf STDERR "time: %s calls: $o->{calls}{$func} $func\n", $o->sec2time( $o->{time}{$func} );
  }
}

sub checkCompleteLinkage {
  my $o = shift;
  my $new_gene = shift;
  my $oid = shift;

  my $gene;
  # compare one gene against all genes in cluster_i (complete linkage)
  # ref $o->{CLUSTERS}[$oid] eq "HASH" or die "\$o->{CLUSTERS}[$oid]=$o->{CLUSTERS}[$oid] not hash ref\n";
  # ref $o->{CLUSTERS}[$oid]{links} eq "HASH" or die "\$o->{CLUSTERS}[$oid]{links}=$o->{CLUSTERS}[$oid]{links} not hash ref\n";
  foreach $gene ( keys %{$o->{CLUSTERS}[$oid]{links}} ) {
    exists $o->{BEST_HIT}{$gene}{$new_gene} and exists $o->{BEST_HIT}{$new_gene}{$gene} and next;
    $o->{VERBOSE} and printf STDERR "INFO: checkCompleteLinkage $new_gene is reciprical with some but not all genes in cluster=$oid (%s)\n", join " ", keys %{$o->{CLUSTERS}[$oid]{links}};
    return;
  }

  $o->{VERBOSE} and printf STDERR "INFO: checkCompleteLinkage $new_gene added to oid=%d\n", $oid;
  $o->{ID2CLSTR}{$new_gene} = $oid;
  foreach $gene ( keys %{$o->{CLUSTERS}[$oid]{links}} ) {
    $o->{CLUSTERS}[$oid]{links}{$new_gene}{$gene} = undef;
    $o->{CLUSTERS}[$oid]{links}{$gene}{$new_gene} = undef;
  }
}

# see if two clusters can be merged
sub checkCompleteMerge {
  my $o = shift;
  my $cluster_i = shift;
  my $cluster_j = shift;

  # if( exists $o->{nomerge}{$cluster_i}{$cluster_j} ) {
  #   # $o->{VERBOSE} and print STDERR "INFO: checkCompleteMerge oids $cluster_i and $cluster_j failed previous merge.\n";
  #   return undef;
  # }

  # list of current genes from two clusters; sort is for debugging only
  my @genes_i = keys %{$o->{CLUSTERS}[$cluster_i]{links}};
  my @genes_j = keys %{$o->{CLUSTERS}[$cluster_j]{links}};
  my ($gene_x, $gene_y, $gene, $oid, $i, $j);
  # compare all genes between two clusters
  foreach $gene_x ( @genes_i ) {
    foreach $gene_y ( @genes_j ) {
      exists $o->{BEST_HIT}{$gene_x}{$gene_y} and exists $o->{BEST_HIT}{$gene_y}{$gene_x} and next; # reciprocal best hits
      $o->{VERBOSE} and printf STDERR "INFO: checkCompleteMerge oid=%d and oid=%d not merged. $gene_x and $gene_y not reciprical\n", $cluster_i, $cluster_j;
      $o->{nomerge}{$cluster_i}{$cluster_j} = undef;
      $o->{nomerge}{$cluster_j}{$cluster_i} = undef;
      return undef; # $gene_x and $gene_y not reciprical best hits; nothing changes
    }
  }

  $o->{VERBOSE} and printf STDERR "INFO: checkCompleteMerge oid=$cluster_j (@genes_j) merged to oid=$cluster_i (@genes_i)\n";

  foreach $gene_y ( @genes_j ) {
    foreach $gene_x ( @genes_i ) {
      # $o->{CLUSTERS}[$cluster_i]{links}{$gene_x}{$gene_y} = $o->{BEST_HIT}{$gene_x}{$gene_y}[0];
      # $o->{CLUSTERS}[$cluster_i]{links}{$gene_y}{$gene_x} = $o->{BEST_HIT}{$gene_y}{$gene_x}[0];
      $o->{CLUSTERS}[$cluster_i]{links}{$gene_x}{$gene_y} = undef;
      $o->{CLUSTERS}[$cluster_i]{links}{$gene_y}{$gene_x} = undef;
    }
    $o->{ID2CLSTR}{$gene_y} = $cluster_i;
  }
  $o->{CLUSTERS}[$cluster_j] = undef; # cluster_j moved to cluster_i above, undefine it
  # Above we linked genes from two clusters.  Here we are adding the internal links from cluster_j
  for( $i = 0; $i < $#genes_j; $i++ ) {
    for( $j = $i + 1; $j < @genes_j; $j++ ) {
      $o->{CLUSTERS}[$cluster_i]{links}{$genes_j[$i]}{$genes_j[$j]} = undef;
      $o->{CLUSTERS}[$cluster_i]{links}{$genes_j[$j]}{$genes_j[$i]} = undef;
    }
  }
  return $cluster_j;
}

sub printClusters {
  my $o = shift;
  my ($i, $gene, $link, @links, $oid, $bits, $pcnt, $alen );
  @{$o->{CLUSTERS}} = sort { scalar keys %{$b->{links}} <=> scalar keys %{$a->{links}} ||
                       join("",sort keys %{$a->{links}}) cmp join("",sort keys %{$b->{links}}) } @{$o->{CLUSTERS}};

  for ( $i = 0; $i < @{$o->{CLUSTERS}}; $i++ ) {
    my ( @genes, @pcnt_total, @bits_total );
    printf ">%d size=%d", $i + 1, scalar keys %{$o->{CLUSTERS}[$i]{links}};
    # for each gene print gene:median_pcnt:median_bits:median_alen:gene_len
    foreach $gene ( sort keys %{$o->{CLUSTERS}[$i]{links}} ) {
      my ( @bits, @pcnt );
      # The following loops over each match between $gene and other proteins in cluster
      foreach $link ( keys %{$o->{CLUSTERS}[$i]{links}{$gene}} ) {
        push @bits, $o->{BEST_HIT}{$gene}{$link}[0];
        push @pcnt, $o->{BEST_HIT}{$gene}{$link}[1];
      }
      push @bits_total, @bits;
      push @pcnt_total, @pcnt;
      push @genes, sprintf "$gene:%.3g:%.5g:$o->{LEN}{$gene}", median( \@pcnt ), median( \@bits );
    }
    printf " bits_median=%.5g pcnt_median=%.4g %s\n", median( \@bits_total ), median( \@pcnt_total ), join( " ", @genes );
  }

  # Singlets
  $o->{SINGLETS} or next;
  $i = scalar( @{$o->{CLUSTERS}} ) + 1;
  print "# singlets\n";
  foreach $gene ( sort { $o->{LEN}{$b} <=> $o->{LEN}{$a} || $a cmp $b } keys %{$o->{LEN}} ) {
    exists $o->{ID2CLSTR}{$gene} and next; # already in ortholog group
    printf ">$i size=1 $gene:0:0:$o->{LEN}{$gene}\n";
    $i++;
  }
}

sub getSeqLengths {
  my $o = shift;
  my $file = shift;
  my $fh = getFH( $file );
  while( <$fh> ) {
    /^(\d+)\s+(\S+)/ or die "error: in gene length file=$file line=$_\n";
    $o->{LEN}{$2} = $1;
  }
}

# This groups paralogs in each genome by complete linkage.  Each a member of a group is reciprical best hit a gene in another genome then all the paralogs are linked to that gene via:
# $o->{BEST_HIT}{$par[i]}{$other_gene} = [bits, pcnt]
# $o->{BEST_HIT}{$other_gene}{$par[i]} = [bits, pcnt]
# Members of paralog groups are also linked via
# $o->{BEST_HIT}{$par[i]}{$par[j]} = [bits, pcnt]
# $o->{BEST_HIT}{$par[j]}{$par[i]} = [bits, pcnt]

# The main difference between this and the previous groupParalogs_old +
# findOrthosParalogs is that this tries to account for the possibility of all
# paralogs in a group being best hit to a particular other genome.
# Adding paralog links to $o->{BEST_HIT} increases memory usage by about 1/3
# $o->{PARALOGS}{$par1}{$par2} = [bits, pcnt]
# $o->{BEST_HIT}{$par1}{$par2} = [bits, pcnt]
sub groupParalogs {
  my $o = shift;
  my ( $bits, $ok, $other_gene, @genes_in_other_genomes, %genes_in_other_genomes );
  my ( $par1, $par2, %par, @par, @ref, %parpar, $porc, %max_bits, @other_genes, $i );
  # This first block does basically what groupParalogs_old does. It defines
  # paralogs within a genome that are closer to each other than to any other
  # genome.
  my $start = time;
  # If these keys aren't sorted then I believe the results are not reproducible.
  # The following block populate $par{$par1} = [par2, par3, ...] which are lists
  # of paralogs that are more similar to par1 than to any other genome
  foreach $par1 ( sort keys %{$o->{PARALOGS}} ) {
    # If it has no hit to another genome then skip it
    exists $o->{BEST_HIT}{$par1} or next;
    # Find paralogs of $par1 that have no better match to other genome
    foreach $par2 ( sort keys %{$o->{PARALOGS}{$par1}} ) {
      if( not exists $o->{PARALOGS}{$par2}{$par1} ) {
        # Low similarity edge cases, where par1 vs par2 is above threshold, but par2 vs par1 isn't
        delete $o->{PARALOGS}{$par1}{$par2};
        next;
      }
      $bits = $o->{PARALOGS}{$par1}{$par2}[0];
      # $max_bits{$par2} represents the maximum bits found in $o->{BEST_HIT}{$par2}{$other_gene}
      # It will be calculated the first time $par2 is tested
      if( exists $max_bits{$par2} ) {
        if($bits >= $max_bits{$par2}) {
          # par1 and par2 more similar to each other than to anything else
          push @{$par{$par1}}, $par2;
        }
        next; # we don't need to go through $o->{BEST_HIT}{$par2}{$other_gene} for this one.
      } elsif ( not exists $o->{BEST_HIT}{$par2} ) {
        # $par2 has no hits to other genomes
        # In this case we don't want to group it as an ortholog with other genomes.
        # Setting it to a very high value so that it never groups with another paralog.
        $max_bits{$par2} = 1e99;
        next;
      } else {
        @other_genes = keys %{$o->{BEST_HIT}{$par2}};
        $max_bits{$par2} = $o->{BEST_HIT}{$par2}{$other_genes[0]}[0];
        for( $i = 1; $i < @other_genes; $i++ ) {
          if($o->{BEST_HIT}{$par2}{$other_genes[$i]}[0] > $max_bits{$par2}) {
            $max_bits{$par2} = $o->{BEST_HIT}{$par2}{$other_genes[$i]}[0];
          }
        }
      }
      if($bits >= $max_bits{$par2}) {
        # par1 and par2 more similar to each other than to anything else
        push @{$par{$par1}}, $par2;
      }
    }
  }
  if( $o->{DEBUG} ) {
    my $fh = getFH( "pars_hash.txt", { mode => ">" } );
    foreach $par1 ( sort keys %par ) {
      printf $fh "$par1 %s\n", join " ", sort @{$par{$par1}};
    }
  }
  $o->logDurationMemory( "groupParalogs get unique paralogs", $o->sec2time(time - $start) ); $start = time;
  # Now cluster these "true" paralogs into groups within each genome; sort of
  # like ortholog clustering. This should save time.
  my $par_clusters = $o->getParalogClusters( \%par );
  $o->logDurationMemory( "groupParalogs getParalogClusters", $o->sec2time(time - $start) ); $start = time;
  # paralog_clusters is for debugging
  # my $fh_porc = getFH( "paralog_clusters", { mode => ">" } );
  # The sort here is for debugging determinism
  foreach $porc ( sort { scalar keys %{$b->{links}} <=> scalar keys %{$a->{links}} } @$par_clusters ) {
    # get best hits for $par1 and @{$par{$par1}}
    # sort is for debugging determinism.  It will affect what @ref is used below
    # @par is list of paralogs in group $porc
    @par = sort keys %{$porc->{links}};
    # printf $fh_porc "%d %s\n", $#par + 1, join " ", @par;
    # @genes_in_other_genomes is a list of best hits for all genes in @par
    @genes_in_other_genomes = ();
    foreach $par1 ( @par ) {
      push @genes_in_other_genomes, keys %{$o->{BEST_HIT}{$par1}};
    }
    # %genes_in_other_genomes is a list of genes from other genomes that are best hits to one or more of the paralogs
    %genes_in_other_genomes = ();
    @genes_in_other_genomes{@genes_in_other_genomes} = ();
    @genes_in_other_genomes = sort keys %genes_in_other_genomes;
    # printf $fh_porc "genes_in_other_genomes: %s\n", join " ", @genes_in_other_genomes;
    foreach $other_gene ( @genes_in_other_genomes ) {
      $ok = 0;
      foreach $par1 ( @par ) {
        # See if one of the paralogs is the best hit for some other genome_gene.
        # If it matches one make it match all of them
        if( exists $o->{BEST_HIT}{$other_gene}{$par1} ) {
          # @ref used as proxy for others since we haven't saved all scores between {$other_gene}{$x}
          @ref = @{$o->{BEST_HIT}{$other_gene}{$par1}}; # (bits, pcnt)
          $ok = 1; # one of that paralogs matched $other_gene
          last;
        }
      }
      if( not $ok ) { # we  didn't find a reciprical link between $other_gene and these paralogs
        next;
      }
      # This is where we actually associate paralogs with other genomes. At
      # least one of the paralogs was best hit for $other_gene so we say all the
      # paralogs are, essentially treating all paralogs as a single gene.
      foreach $par1 ( @par ) {
        if( not exists $o->{BEST_HIT}{$other_gene}{$par1} ) {
          $o->{BEST_HIT}{$other_gene}{$par1} = [ @ref ];
        }
        if( not exists $o->{BEST_HIT}{$par1}{$other_gene} ) {
          $o->{BEST_HIT}{$par1}{$other_gene} = [ @ref ];
        }
        # If $other_gene also has paralogs then the reciprical populating will
        # be done on {BEST_HIT}{$par1}{$other_gene} etc
      }
    }
    # Add links between paralogs
    foreach $par1 ( @par ) {
      foreach $par2 ( @par ) {
        $par2 eq $par1 and next;
        $o->{BEST_HIT}{$par1}{$par2} = [ @{$o->{PARALOGS}{$par1}{$par2}} ];
      }
    }
  } # end of foreach $porc ( @$par_clusters )
  $o->logDurationMemory( "groupParalogs over \@\$par_clusters", $o->sec2time(time - $start) ); $start = time;
  delete $o->{PARALOGS};
}

sub getParalogClusters {
  my $o = shift;
  my $pars = shift;
  my ( %hits, $q_par, @pars, $i, $j, %porthos, $G );
  # The "sorts" here is to make it
  # deterministic. $o->{PARALOGS}{$pars[$i]}{$pars[$j]} can have slightly
  # different scores than the reciprical which can make the median pcnt and bit
  # values slightly different in the final output and can also slightly alter
  # final groupings. Making it deterministic is mainly for debugging; to be able
  # to more easily see what differences are due to changes in the code and not
  # randomness due to random hash key lists.

  # This block populates $porthos{$G} where $G represents a particular
  # genome. For each $G porthos is a list of paralog pairs and their associated
  # similarity scores. porthos is used to cluster paralogs into groups within
  # each genome.
  foreach $q_par ( sort keys %$pars ) {
    ( $G = $q_par ) =~ s/_.*//;
    @pars = sort ( $q_par, @{$pars->{$q_par}} );
    for( $i = 0; $i < $#pars; $i++ ) {
      for( $j = $i + 1; $j < @pars; $j++ ) {
        # paralogs with genome $G
        exists $hits{$pars[$i]}{$pars[$j]} and next;
        # exists $hits{$pars[$j]}{$pars[$i]} and next; # Is this needed?
        # $done{$pars[$i]}{$pars[$j]} = undef;
        # $done{$pars[$j]}{$pars[$i]} = undef;
        # The following is slightly hard to understand, but it's possible
        # $o->{PARALOGS}{$pars[$i]}{$pars[$j]} doesn't exist if they are members
        # of @pars with a different query like $pars->{par_k} = [par_i,
        # par_j]. par_k could be above threshold to both par_i and par_j, but
        # par_i and par_j may be below threshold to each other.
        exists $o->{PARALOGS}{$pars[$i]}{$pars[$j]} or next;
        # exists $o->{PARALOGS}{$pars[$j]}{$pars[$i]} or next;
        # This puts all paralog pairs into same taxon array.
        # They are clustered below.
        push @{$porthos{$G}}, [ $o->{PARALOGS}{$pars[$i]}{$pars[$j]}[0], $pars[$i], $pars[$j] ];
        $hits{$pars[$i]}{$pars[$j]} = undef;
        $hits{$pars[$j]}{$pars[$i]} = undef;
      }
    }
  }
  my ( $por, $bits, $idA, $idB, %id2clstr, %nomerge, @cluster, @ga, @gb, $removed, @removed );
  my $cnt = 1;
  # my $ofh = getFH( "porthos.txt", { mode => ">" } );
  foreach $G ( sort { @{$porthos{$b}} <=> @{$porthos{$a}} || $a <=> $b } keys %porthos ) {
    foreach $por ( sort { $b->[0] <=> $a->[0] } @{$porthos{$G}} ) { # high bits first
      ($bits, $idA, $idB) = @$por;
      # print $ofh "$G $bits $idA $idB\n";
      if ( exists $id2clstr{$idA} and exists $id2clstr{$idB} ) {
        # both already in a cluster
        exists $nomerge{$id2clstr{$idA}}{$id2clstr{$idB}} and next;
        $removed = checkCompleteMerge_paralog( $id2clstr{$idA}, $id2clstr{$idB}, \%id2clstr, \@cluster, \%hits );
        if( defined $removed ) {
          # $removed is $id2clstr{$idB} which was merged with $id2clstr{$idA}
          push @removed, $removed; # for cleaning the @cluster at the end
        } else {
          # mark these clusters as unable to be merged so they don't need to be checked again
          $nomerge{$id2clstr{$idA}}{$id2clstr{$idB}} = undef;
          $nomerge{$id2clstr{$idB}}{$id2clstr{$idA}} = undef;
        }
      } elsif( exists $id2clstr{$idA} ) {
        # Try to add $idB to 
        checkCompleteLinkage_paralog( $idB, $id2clstr{$idA}, $cluster[$id2clstr{$idA}], \%hits, \%id2clstr );
      } elsif( exists $id2clstr{$idB} ) {
        checkCompleteLinkage_paralog( $idA, $id2clstr{$idB}, $cluster[$id2clstr{$idB}], \%hits, \%id2clstr );
      } else {
        # create a new cluster
        $id2clstr{$idA} = $cnt;
        $id2clstr{$idB} = $cnt;
        $cluster[$cnt]{links}{$idA}{$idB} = undef;
        $cluster[$cnt]{links}{$idB}{$idA} = undef;
        $cnt++;
        # while( defined $cluster[$cnt] ) { $cnt++ }
      }
    }
  }
  foreach $cnt ( sort { $b <=> $a } @removed ) {
    defined $cluster[$cnt] or splice @cluster, $cnt, 1;
    print STDERR "DEBUG: \@cluster[$cnt] spliced\n";
  }
  return \@cluster;
}

sub checkCompleteLinkage_paralog {
  my $new_gene = shift;
  my $oid = shift;
  my $cluster = shift;
  my $hits = shift;
  my $id2clstr = shift;
  my ( $gene );
  foreach $gene ( keys %{$cluster->{links}} ) {
    exists $hits->{$gene}{$new_gene} and exists $hits->{$new_gene}{$gene} and next;
    return undef; # something didn't link; no merge
  }
  $id2clstr->{$new_gene} = $oid;
  foreach $gene ( keys %{$cluster->{links}} ) {
    $cluster->{links}{$new_gene}{$gene} = undef;
    $cluster->{links}{$gene}{$new_gene} = undef;
  }
}

sub checkCompleteMerge_paralog {
  my $ca = shift;
  my $cb = shift;
  my $id2clstr = shift;
  my $cluster = shift;
  my $hits = shift;
  my @ga = keys %{$cluster->[$ca]};
  my @gb = keys %{$cluster->[$cb]};
  my ( $ga, $gb, $i, $j );
  foreach $ga ( @ga ) {
    foreach $gb ( @gb ) {
      # This checks that all genes in clusters a and b are linked
      exists $hits->{$ga}{$gb} and exists $hits->{$gb}{$ga} and next;
      return undef;
    }
  }
  # Merge cluster b into cluster a
  foreach $gb ( @gb ) {
    foreach $ga ( @ga ) {
      $cluster->[$ca]{links}{$ga}{$gb} = undef;
      $cluster->[$ca]{links}{$gb}{$ga} = undef;
    }
    $id2clstr->{$gb} = $ca;
  }
  $cluster->[$cb] = undef;
  # add cb vs cb links
  for( $i = 0; $i < $#gb; $i++ ) {
    for( $j = $i + 1; $j < @gb; $j++ ) {
      $cluster->[$ca]{links}{$gb[$i]}{$gb[$j]} = undef;
      $cluster->[$ca]{links}{$gb[$j]}{$gb[$i]} = undef;
    }
  }
  return $cb;
}

# This was part of the old way in combination with findOrthosParalogs.  This is
# faster than than the current groupParalogs but may not be as accurate.
sub groupParalogs_old {
  my $o = shift;
  my ( $qid, $hid, $bits, $ok, $bh );
  foreach $qid ( sort keys %{$o->{PARALOGS}} ) {
    foreach $hid ( sort keys %{$o->{PARALOGS}{$qid}} ) {
      $bits = $o->{PARALOGS}{$qid}{$hid}[0];
      $ok = 1;
      foreach $bh ( keys %{$o->{BEST_HIT}{$hid}} ) { # is $hid more similar to a gene in some other genome?
        if( $o->{BEST_HIT}{$hid}{$bh}[0] > $bits ) {
          $ok = 0;
          last;
        }
      }
      $ok or next;
      # $hid has best hit to $qid relative to other genomes
      if( exists $o->{PARALOGS}{$hid}{$qid} ) { # reciprical paralog
        @{$o->{BEST_HIT}{$hid}{$qid}} = @{$o->{PARALOGS}{$hid}{$qid}};
        push @{$o->{QID2PARALOGS}{$qid}}, $hid;
      }
    }
  }
  delete $o->{PARALOGS};
}

# Populating @{$o->{ORTHOS}} increases memory usage by about 1/3
sub findOrthos {
  my $o = shift;
  my ($qid, $hid, %done);
  foreach $qid ( sort keys %{$o->{BEST_HIT}} ) {
    foreach $hid ( sort keys %{$o->{BEST_HIT}{$qid}} ) {
      if( not exists $o->{BEST_HIT}{$hid}{$qid} ) {
        delete $o->{BEST_HIT}{$qid}{$hid}; # not reciprical
        next;
      }
      exists $done{$qid}{$hid} and next;
      $done{$hid}{$qid} = (); # so we don't duplicate ortho in @{$o->{ORTHOS}} array
      push @{$o->{ORTHOS}}, [ $o->{BEST_HIT}{$qid}{$hid}[0], $qid, $hid ]; # [score, gene_i, gene_j]
    }
  }
  # highest scoring pairs first
  @{$o->{ORTHOS}} = sort { $b->[0] <=> $a->[0] || $a->[1] cmp $b->[1] || $a->[2] cmp $b->[2] } @{$o->{ORTHOS}};
}

# Write all the reciprical best hits
sub writeRecipricalBestHits {
  my $o = shift;
  my ($qid, $hid);
  my $fh_hits = getFH( "rbh.dat", { mode => ">" } );
  foreach $qid ( sort keys %{$o->{BEST_HIT}} ) {
    foreach $hid ( sort keys %{$o->{BEST_HIT}{$qid}} ) {
      printf $fh_hits "$qid\t$hid\t%s\n", join "\t", @{$o->{BEST_HIT}{$qid}{$hid}};
    }
  }
}

# Deprecated. Superseded by groupParalogs
sub findOrthosParalogs {
  my $o = shift;
  my ( $qid, $hid, %done, $paralog, $rep );
  foreach $qid (sort keys %{$o->{BEST_HIT}} ) {
    $done{$qid}{$qid} = (); # To prevent this due to $qid -> $hid -> $qid link.
    foreach $hid (sort keys %{$o->{BEST_HIT}{$qid}}) {
      if( not exists $o->{BEST_HIT}{$hid}{$qid} ) { # not reciprical, but maybe paralog if $hid is
        $rep = undef;
        if( exists $o->{QID2PARALOGS}{$hid} ) {
          foreach $paralog ( @{$o->{QID2PARALOGS}{$hid}} ) {
            if( exists $o->{BEST_HIT}{$paralog}{$qid} ) {
              # If one matches they all match because they are defined as paralogs that are specific to that genome.
              # It is conceivable that some genes in this group, more distant to $paralog should go into different homologous
              # ortholog group, but for now we are just putting them all in the same ortholog group.
              $rep = $paralog;
              last;
            }
          }
        }
      } else { # normal reciprical best hit
        $rep = $hid;
      }
      defined $rep or next; # $hid not reciprical and no paralogs match
      foreach $paralog ( @{$o->{QID2PARALOGS}{$hid}} ) {
        exists $done{$qid}{$paralog} and next;
        exists $done{$paralog}{$qid} and next;
        $done{$paralog}{$qid} = ();
        # since $qid/$paralog could come up again from other paralogs or as $hid we need to mark it as done
        $done{$qid}{$paralog} = ();
        # For printing clusters we need to put something in %o->{BEST_HIT} if it isn't there
        # Just inheriting from representative paralog
        # make $qid reciprical with all paralogs
        if( not exists $o->{BEST_HIT}{$paralog}{$qid} ) {
          $o->{BEST_HIT}{$paralog}{$qid} = [ @{$o->{BEST_HIT}{$rep}{$qid}} ];
        }
        $o->{BEST_HIT}{$qid}{$paralog} = [ @{$o->{BEST_HIT}{$qid}{$hid}} ];
        # Giving all paralogs $o->{BEST_HIT}{$qid}{$rep}{bits}. Possibly better choice? Probably doesn't matter?
        $qid eq $paralog and die "error: same gene qid=$qid rep=$rep\n";
        push @{$o->{ORTHOS}}, [ $o->{BEST_HIT}{$qid}{$hid}[0], $qid, $paralog];
      }
      exists $done{$qid}{$hid} and next;
      exists $done{$hid}{$qid} and next;
      $done{$hid}{$qid} = ();
      $done{$qid}{$hid} = ();
      if( not exists $o->{BEST_HIT}{$hid}{$qid} ) {
        $o->{BEST_HIT}{$hid}{$qid} = [ @{$o->{BEST_HIT}{$rep}{$qid}} ];
      }
      $qid eq $hid and die "error: same gene qid=$qid rep=$rep hid=$hid 2\n";
      push @{$o->{ORTHOS}}, [ $o->{BEST_HIT}{$qid}{$hid}[0], $qid, $hid];
    }
  }
  @{$o->{ORTHOS}} = sort {$b->[0] <=> $a->[0]} @{$o->{ORTHOS}}; # highest scoring pairs first
}

# This is mainly for estimating overall bit score and number of identities
# between a query and subjct that is composed of potentially overlapping hsps.


# HSP filter.  If the local percent identity of an HSP differs from the first
# HSP by more than $diff_pcnt then it is skipped.
# $hr can currently contain $hr->{diff_pcnt} but more options can be added.
# diff_pcnt should be between 0 and 1

# $hsp = [ $bits, $nident, $qb, $qe, $sb, $se ];
sub mergeHsps {
  my $hsps = shift; # hsps for a given query/sbjct pair
  my $hr = shift;
  my ( $qcov, $scov, $qcovNew, $scovNew, $newFrac, $i, @qcov, @scov );
  defined $hr or $hr = {};
  my $diff_pcnt = defined $hr->{diff_pcnt} ? $hr->{diff_pcnt} : 1e9;
  @$hsps = sort { $b->[0] <=> $a->[0] } @$hsps; # highest bits first
  my $hsp = $hsps->[0]; # shift @$hsps;
  my $bits = $hsp->[0];
  my $nident = $hsp->[1];
  my $nhsps = 1;
  my $qcovTotal = $hsp->[3] - $hsp->[2] + 1;
  my $scovTotal = $hsp->[5] - $hsp->[4] + 1;
  my $pcntFirst = $nident / $qcovTotal;
  @qcov[ $hsp->[2]..$hsp->[3] ] = (1) x $qcovTotal;
  @scov[ $hsp->[4]..$hsp->[5] ] = (1) x $scovTotal;
  foreach $hsp ( @{$hsps}[1..$#$hsps] ) {
    $qcov = $hsp->[3] - $hsp->[2] + 1;
    $scov = $hsp->[5] - $hsp->[4] + 1;
    # Using absolute value here because we don't want things with much lower
    # identity to be included, but we also don't want some short high identity
    # hsp to be combined with a longer lower identity region as they are
    # unlikely to correspond to the same alignment.
    if( abs( $hsp->[1] / $qcov - $pcntFirst ) > $diff_pcnt ) {
      next;
    }
    ($qcovNew, $scovNew) = (0, 0);
    for($i = $hsp->[2]; $i <= $hsp->[3]; $i++) {
      defined $qcov[$i] and next;
      $qcov[$i] = 1;
      $qcovNew++;
    }
    for($i = $hsp->[4]; $i <= $hsp->[5]; $i++) {
      defined $scov[$i] and next;
      $scov[$i] = 1;
      $scovNew++;
    }
    # The min() here is because if either the query or the sbjct is reused we don't want to include it.
    $newFrac = min($qcovNew / $qcov, $scovNew / $scov ); # fraction of HSP that is unique coverage
    $bits += $hsp->[0] * $newFrac;
    $nident += $hsp->[1] * $newFrac;
    $qcovTotal += $qcovNew;
    $scovTotal += $scovNew;
    $nhsps++;
  }
  return ( $bits, $nident, $qcovTotal, $scovTotal );
}

# Can handle an array of numbers and array_refs mixed, e.g. max( 8, [5, 6, 7], 99, [13] )
# Be sure not to modify array references
sub max {
  @_ or croak "error: DnaNinja::Utils::max - no arguments provided";
  # Initialize $max
  my ( $max, $val, $refOrVal );
  if ( ref $_[0] ) {
    defined $_[0][0] or croak "error: DnaNinja::Utils::max - first argument passed was an empty array reference.";
    $max = $_[0][0];
  } else {
    $max = shift;
  }

  foreach $refOrVal ( @_ ) {
    if ( ref $refOrVal ) {
      foreach $val ( @$refOrVal ) {
        $val > $max and $max = $val;
      }
    } else {
      $refOrVal > $max and $max = $refOrVal;
    }
  }
  return $max;
}

sub min {
  my $min = shift(@_);
  foreach my $foo (@_) {
    $foo < $min and $min = $foo;
  }
  return $min;
}

# WARNING: This used to check all files for reading with `file -L '$file'`. This
# was very slow so I have now removed it. Now if you are reading a file that is
# compressed you can specify the program to uncompress it with the attribute
# {unzip => program} where "program" can be "zcat", "bzcat", etc. Otherwise it
# will uncompress a file ending in .gz with zcat.
sub getFH {
  my $file = shift;
  my $hr = shift;
  defined $file or croak "error: DnaNinja::Utils::getFH - file variable not defined";
  my %opts = (mode => "<", seekable => 0, maxsize => 5e9, unzip => undef );
  if( defined $hr ) {
    ref $hr eq "HASH" or croak "error: DnaNinja::Utils::getFH second argument of getFH must now be a hash ref with possible keys {mode, seekable, and maxsize}";
    my $key;
    foreach $key( qw(mode seekable maxsize unzip) ) {
      exists $hr->{$key} or next;
      $opts{$key} = $hr->{$key};
      delete $hr->{$key};
    }
    if( %$hr ) {
      print STDERR "ERROR: unrecognized keys in hash reference passed to getFH\n";
      foreach $key ( keys %$hr ) {
        print STDERR "$key = $hr->{$key}\n";
      }
      exit 1;
    }
  }

  my ( $fh );
  if ( $file eq '-' ) {
    $opts{mode} =~ />/ and return \*STDOUT;
    $opts{seekable} or return \*STDIN;
    return seekable( \*STDIN );
  }
  if( $opts{mode} =~ />/ ) {
    open $fh, $opts{mode}, $file or die "error: DnaNinja::Utils::getFH can't open file=$file in mode '$opts{mode}'\n";
    return $fh;
  }

  # defined $opts{unzip} or $opts{unzip} = getUnzip( $file );
  my $unzip;
  if ( defined $opts{unzip} ) {
    open $fh, "$opts{unzip} '$file'|" or croak "ERROR: DnaNinja::Utils::getFH can't open $opts{unzip} '$file'|";
    $opts{seekable} and return seekable( $fh );
    return $fh;
    # Disabling this because it is very slow
    # } elsif ($unzip = getUnzip($file)) {
  } elsif ($file =~ /\.gz$/) {
    open $fh, "zcat '$file'|" or croak "ERROR: DnaNinja::Utils::getFH can't open zcat '$file'|";
    $opts{seekable} and return seekable( $fh );
    return $fh;
  } else {
    open $fh, $opts{mode}, "$file" or croak "getFH: DnaNinja::Utils::getFH error: can't open file=$file in mode=$opts{mode}";
    return $fh;
  }
  croak "ERROR: DnaNinja::Utils::getFH( $file ) with mode => $opts{mode}, seekable => $opts{seekable}, maxsize => $opts{maxsize}";
}

# NOTE: I am no longer using an in-memory file. Doing regular expression matches
# on lines from in-memory files is very very slow on current cygwin build so avoid using them.
# I posted a bug for this at https://github.com/Perl/perl5/issues/21877
sub seekable {
  my $fh = shift;
  my $fh2 = tempfile(UNLINK => 1);
  while( <$fh> ) {
    print $fh2;
  }
  seek($fh2, 0, 0);
  return $fh2;
}

sub mean {
  my $ar;
  if ( @_ > 1 ) {
    $ar = [@_];
  } else {
    $ar = shift;
  }
  my $sum = 0;
  for (my $i = 0; $i < @$ar; $i++) {
    $sum += $ar->[$i]
  }
  return $sum / @$ar;
}

# can pass a hash ref with key=even_rule to determine how to define median for an even number of data
# possible values are low, high and mean.  Default is mean.
sub median {
  my $ar = shift;

  # If odd, no need to check even_rule
  my @s = sort {$a <=> $b} @$ar;
  my $n = @s;
  if ($n % 2) {
    return $s[int $n / 2];
  }

  my $h = shift;
  my $even_rule = 'mean';
  if (exists $h->{even_rule}) {
    $even_rule = "\L$h->{even_rule}";
    $even_rule =~ /^low|high|mean$/ or die "error: Statistics::median even_rule must be low high or mean\n";
  }
  if ($even_rule eq 'mean') {
    return ($s[$n / 2 - 1] + $s[$n / 2]) / 2;
  } elsif ($even_rule eq 'high') {
    return $s[$n / 2];
  } else {  # low
    return $s[$n / 2 - 1];
  }
}

1;
