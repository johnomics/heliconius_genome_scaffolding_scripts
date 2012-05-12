#!/usr/bin/env perl

# scaffold_heliconius_genome.pl

# This script is provided FOR INFORMATION ONLY and is not intended for use.
# It was written specifically for the Heliconius melpomene genome paper
# (Heliconius Genome Consortium, doi: 10.1038/nature11041) and has not
# been adapted for general use.

# Purpose: Synthesise linkage map of H. melpomene with genome scaffolds
#          to produce a chromosomal AGP file,
#          revising scaffolds where they feature on multiple linkage groups

# Input  : marker position list from convert_vcf_to_joinmap_heliconius.pl
#          lmxll markers (chromosome prints)
#          scaffold AGP file
#          gene annotation GFF file
#          folder of linkage maps produced from markers output by
#            convert_vcf_to_joinmap_heliconius.pl
#          gene marker file listing scaffolds for genes with known
#            chromosomal positions, for consistency with previous maps
#            eg a known gene location is listed in CSV format:
#              1,Beta-actin,scf0000001
#            this will link the linkage group containing scf0000001
#            to chromosome 1

# Output : revised scaffold AGP file
#          chromosome AGP file
#          full scaffolding details and statistics

# Author: John Davey john.davey@ed.ac.uk
# Begun 30/08/11

#############################################################################
###                                                                       ###
###                                 CODE                                  ###
###                                                                       ###
#############################################################################

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Data::Dumper;
use List::Util qw/max/;

# Autoflush output so reporting on progress works
$| = 1;

my $marker_filename = "";
my $female_filename = "";
my $lg_filename = "";
my $agp_filename = "";
my $gff_filename = "";
my $maps_folder  = "";
my $scaffold_agp_output   = "";
my $chromosome_agp_output = "";
my $debug;
my $verbose;

my $options_okay = GetOptions(
    'markers=s'    => \$marker_filename,
    'female=s'     => \$female_filename,
    'lg=s'         => \$lg_filename,
    'agp=s'        => \$agp_filename,
    'gff=s'        => \$gff_filename,
    'joinmap=s'    => \$maps_folder,
    'scaffold=s'   => \$scaffold_agp_output,
    'chromosome=s' => \$chromosome_agp_output,
    'debug'        => \$debug,
    'verbose'      => \$verbose,
);

croak "No markers file! Please specify -m $OS_ERROR\n"
  if ( $marker_filename eq "" );
croak "No female markers file! Please specify -f $OS_ERROR\n"
  if ( $female_filename eq "" );
croak "No gene marker file! Please specify -l $OS_ERROR\n"
  if ( $lg_filename eq "" );
croak "No AGP file! Please specify -a $OS_ERROR\n" if ( $agp_filename eq "" );
croak "No GFF file! Please specify -g $OS_ERROR\n" if ( $gff_filename eq "" );
croak "No maps folder! Please specify -j $OS_ERROR\n" if ( $maps_folder eq "" );

croak "No scaffold AGP output filename given! Please specify -s $OS_ERROR\n"
  if ( $scaffold_agp_output eq "" );
croak "No chromosome AGP output filename given! Please specify -c $OS_ERROR\n"
  if ( $chromosome_agp_output eq "" );

my %scf_contigs;
my %scf_class;
my %scf_lengths;
my %revised_genome;

my $cur_scf           = "";
my $scf_length        = 0;
my $scf_length_no_gap = 0;
open my $agp_file, '<', $agp_filename
  or croak "Can't open $agp_filename $OS_ERROR!\n";
while ( my $agp_line = <$agp_file> ) {
    chomp $agp_line;
    my @agp_fields  = split /\t/, $agp_line;
    my $scf         = $agp_fields[0];
    my $contig_num  = $agp_fields[3];
    my $contig_type = $agp_fields[4];

    if ( $cur_scf ne $scf ) {
        if ( $cur_scf ne "" ) {
            $scf_lengths{gap}{$cur_scf}              = $scf_length;
            $scf_lengths{no_gap}{$cur_scf}           = $scf_length_no_gap;
            $revised_genome{length}{$cur_scf}        = $scf_length;
            $revised_genome{length_no_gap}{$cur_scf} = $scf_length_no_gap;

        }
        $cur_scf = $scf;
        $scf_class{all}{$scf}++;
        $scf_class{revised}{$scf}++;
        $scf_length_no_gap = 0;
    }
    $scf_contigs{$scf}{$contig_num}{scf_start} = $agp_fields[1];
    $scf_contigs{$scf}{$contig_num}{scf_end}   = $agp_fields[2];
    $scf_length = $scf_contigs{$scf}{$contig_num}{scf_end};

    $scf_contigs{$scf}{$contig_num}{type} = $contig_type;
    if ( $contig_type eq "W" ) {
        $scf_length_no_gap +=
          $scf_contigs{$scf}{$contig_num}{scf_end} -
          $scf_contigs{$scf}{$contig_num}{scf_start} + 1;
        $scf_contigs{$scf}{$contig_num}{name}        = $agp_fields[5];
        $scf_contigs{$scf}{$contig_num}{ctg_start}   = $agp_fields[6];
        $scf_contigs{$scf}{$contig_num}{ctg_end}     = $agp_fields[7];
        $scf_contigs{$scf}{$contig_num}{orientation} = $agp_fields[8];
    }
    elsif ( $contig_type eq "N" ) {
        $scf_contigs{$scf}{$contig_num}{length} = $agp_fields[5];
    }
}
close $agp_file;
$scf_lengths{gap}{$cur_scf}              = $scf_length;
$scf_lengths{no_gap}{$cur_scf}           = $scf_length_no_gap;
$revised_genome{length}{$cur_scf}        = $scf_length;
$revised_genome{length_no_gap}{$cur_scf} = $scf_length_no_gap;

my %scf_genes;
open my $gff_file, '<', $gff_filename
  or croak "Can't open $gff_filename $OS_ERROR!\n";
while ( my $gff_line = <$gff_file> ) {
    next if ( $gff_line !~ /gene/ );
    chomp $gff_line;
    my @gff_f     = split /\t/, $gff_line;
    my $scf       = $gff_f[0];
    my $gene_name = $gff_f[8];
    if ( $gene_name =~ /ID=(.+?);/ ) {
        $gene_name = $1;
    }
    $scf_genes{$scf}{$gene_name}{start} = $gff_f[3];
    $scf_genes{$scf}{$gene_name}{end}   = $gff_f[4];
}
close $gff_file;

my %map_markers;
open my $female_file, '<', $female_filename
  or croak "Can't open $female_filename $OS_ERROR!\n";
while ( my $female_line = <$female_file> ) {
    chomp $female_line;
    if ( $female_line =~ /\<lmxll\>/ ) {
        my ( $marker, $type, $genotypes ) = split /\t/, $female_line;
        my ( $id, $pos, $hkhk, $ll_head ) = split /[_:]/, $marker;
        $map_markers{$id}{pos}   = $pos;
        $map_markers{$id}{hkhk}  = $hkhk;
        $map_markers{$id}{chrpr} = $ll_head;
        $map_markers{$id}{gt}    = $genotypes;
    }
}
close $female_file;

my %chrpr_to_lg;
my %map;
my @maps_filenames = glob("$maps_folder/*");
foreach my $maps_filename (@maps_filenames) {
    open my $maps_file, '<', $maps_filename
      or croak "Can't open $maps_filename $OS_ERROR!\n";
    my $group_num = 0;
    my $chrpr     = 0;

    my $cm_order = 0;
    my $cm       = -1;

    my @non_ll_markers;
    while ( my $maps_line = <$maps_file> ) {
        next if ( $maps_line =~ /^;/ );
        next if ( $maps_line eq "" );
        if ( $maps_line =~ /^group (.+)$/ ) {
            $group_num = $1;
        }
        elsif ( $maps_line =~ /^(.+?)\s+(.+)$/ ) {
            my @id_fields = split /[_:]/, $1;
            my $id = $id_fields[0];
            if ( $2 ne $cm ) {
                $cm = $2;
                $cm_order++;
            }
            $map_markers{$id}{pos}      = $id_fields[1];
            $map_markers{$id}{cm}       = $cm;
            $map_markers{$id}{cm_order} = $cm_order;

            $map{$group_num}{$cm}{$id}++;

            if ( ( defined $id_fields[2] ) && ( $id_fields[2] ne 'sex' ) ) {
                $map_markers{$id}{chrpr} = $id_fields[2];
                $chrpr = $id_fields[2];
            }
            else {
                push @non_ll_markers, $id;
            }
        }
    }
    close $maps_file;
    $chrpr_to_lg{$chrpr} = $group_num;
    map { $map_markers{$_}{chrpr} = $chrpr; } @non_ll_markers;
}

# Load scaffold positions using revised lmxll markers
open my $scf_marker_file, '<', $marker_filename
  or croak "Can't open $marker_filename $OS_ERROR!\n";

my %all_scf_marker;
my %map_scf_marker;
while ( my $marker_line = <$scf_marker_file> ) {
    chomp $marker_line;
    my (
        $marker,  $scfpos,  $type,      $genotypes,
        $depth,   $qual,    $pos_dp,    $mq,
        $male_dp, $male_gq, $female_dp, $female_gq
      )
      = split /\t/,
      $marker_line;
    next if ( ( $pos_dp > 1000 ) && ( $mq <= 66.54 ) );
    next if ( $pos_dp >= 3631 );
    my ( $scf, $pos ) = split /:/, $scfpos;
    $all_scf_marker{$scf}++;
    $scf_class{markers}{$scf}++;
    if ( defined $map_markers{$marker} ) {
        $map_scf_marker{$scf}{$marker}{$pos}++;
        $map_markers{$marker}{gt} = $genotypes;

        $scf_class{map}{$scf}++;
        foreach my $contig ( keys %{ $scf_contigs{$scf} } ) {
            if (   ( $scf_contigs{$scf}{$contig}{scf_start} <= $pos )
                && ( $scf_contigs{$scf}{$contig}{scf_end} >= $pos ) )
            {
                $map_markers{$marker}{scf}{$scf}{$contig}++;
                $scf_contigs{$scf}{$contig}{markers}{$pos}{type}    = $type;
                $scf_contigs{$scf}{$contig}{markers}{$pos}{marker}  = $marker;
                $scf_contigs{$scf}{$contig}{markers}{$pos}{depth}   = $depth;
                $scf_contigs{$scf}{$contig}{markers}{$pos}{qual}    = $qual;
                $scf_contigs{$scf}{$contig}{markers}{$pos}{pos_dp}  = $pos_dp;
                $scf_contigs{$scf}{$contig}{markers}{$pos}{mq}      = $mq;
                $scf_contigs{$scf}{$contig}{markers}{$pos}{male_dp} = $male_dp;
                $scf_contigs{$scf}{$contig}{markers}{$pos}{male_gq} = $male_gq;
                $scf_contigs{$scf}{$contig}{markers}{$pos}{female_dp} =
                  $female_dp;
                $scf_contigs{$scf}{$contig}{markers}{$pos}{female_gq} =
                  $female_gq;
            }
        }
    }
}
close $scf_marker_file;

my %revised_agp;
my $rev_scf = 7190000000001;
map {
    if ( !defined $scf_class{markers}{$_} )
    {
        $scf_class{nomarkers}{$_}++;
    }
} keys %{ $scf_lengths{gap} };

# Assign chromosome prints to contigs wherever possible;
# find genes bridging two contigs and assign chromosome prints to contigs
# with a gene link and a neighbouring contig with a print

foreach my $scf ( keys %map_scf_marker ) {
    my %gene_links;
    foreach my $contig ( sort { $a <=> $b } keys %{ $scf_contigs{$scf} } ) {
        next if ( $scf_contigs{$scf}{$contig}{type} eq "N" );
        foreach my $pos (
            sort { $a <=> $b }
            keys %{ $scf_contigs{$scf}{$contig}{markers} }
          )
        {
            my $marker = $scf_contigs{$scf}{$contig}{markers}{$pos}{marker};
            $scf_contigs{$scf}{$contig}{chrpr}{ $map_markers{$marker}{chrpr} } =
              $pos;
        }
        if ( defined $scf_genes{$scf} ) {
            foreach my $gene ( sort keys %{ $scf_genes{$scf} } ) {
                if (
                    (
                        $scf_genes{$scf}{$gene}{start} <
                        $scf_contigs{$scf}{$contig}{scf_end}
                    )
                    && ( $scf_genes{$scf}{$gene}{end} >
                        $scf_contigs{$scf}{$contig}{scf_start} )
                  )
                {
                    $gene_links{$gene}{$contig}++;
                }
            }
        }
    }
    foreach my $gene ( keys %gene_links ) {
        foreach my $contig1 ( keys %{ $gene_links{$gene} } ) {
            foreach my $contig2 ( keys %{ $gene_links{$gene} } ) {
                next if ( $contig1 eq $contig2 );
                $scf_contigs{$scf}{$contig1}{gene_link}{$contig2}{$gene}++;
            }
        }
    }
}

print
"Scaffold          \tContig\tPos\tMarker\tType\tSupport\tDepth\tQual\tPos_DP\tMapQ\tMaleDP\tFemDP\tDPdiff\tMaleGQ\tFemGQ\tGQdiff\tChrom\tGenotypes\n"
  if $debug;
foreach my $scf (
    sort { $scf_lengths{gap}{$a} <=> $scf_lengths{gap}{$b} }
    keys %map_scf_marker
  )
{
    my %scf_pos_markers;
    my $scf_start_contig    = 1;
    my $scf_end_contig      = 0;
    my $scf_cm_end_contig   = 0;
    my $scf_chrpr           = 0;
    my $candidate_scf_chrpr = 0;
    my $marker_output       = "";
    my %split_scf;
    my $prev_cm_order = "-";
    my $prev_cm_chrpr = 0;

    my %contig_cm_leaps;

    foreach my $contig ( sort { $a <=> $b } keys %{ $scf_contigs{$scf} } ) {
        next if ( $scf_contigs{$scf}{$contig}{type} eq "N" );
        my $emitted_this_contig = 0;

        if ( !defined $scf_contigs{$scf}{$contig}{chrpr} ) {
            my %gene_chrpr;
            foreach
              my $gene_link ( keys %{ $scf_contigs{$scf}{$contig}{gene_link} } )
            {
                if (   ( defined $scf_contigs{$scf}{$gene_link}{chrpr} )
                    && ( keys %{ $scf_contigs{$scf}{$gene_link}{chrpr} } == 1 )
                  )
                {
                    my $gene_link_chrpr =
                      ( keys %{ $scf_contigs{$scf}{$gene_link}{chrpr} } )[0];
                    foreach my $gene (
                        keys
                        %{ $scf_contigs{$scf}{$contig}{gene_link}{$gene_link} }
                      )
                    {
                        $gene_chrpr{$gene_link_chrpr}{$gene}++;
                        $gene_chrpr{$gene_link_chrpr}{$gene}++;
                    }
                }
            }
            if ( keys %gene_chrpr == 1 ) {
                my $linked_chrpr = ( keys %gene_chrpr )[0];
                print "$scf\t$contig\t$linked_chrpr" if $debug;

                foreach my $gene ( keys %{ $gene_chrpr{$linked_chrpr} } ) {
                    $scf_contigs{$scf}{$contig}{chrpr}{$linked_chrpr}{$gene}++;
                    print "\t$gene" if $debug;

                }
                print "\n" if $debug;
            }
        }

        my %contig_cm;
        my %contig_cm_order;

        foreach my $pos (
            sort { $a <=> $b }
            keys %{ $scf_contigs{$scf}{$contig}{markers} }
          )
        {
            my $sex_dp_diff =
              abs( $scf_contigs{$scf}{$contig}{markers}{$pos}{female_dp} -
                  $scf_contigs{$scf}{$contig}{markers}{$pos}{male_dp} );
            my $sex_gq_diff =
              abs( $scf_contigs{$scf}{$contig}{markers}{$pos}{female_gq} -
                  $scf_contigs{$scf}{$contig}{markers}{$pos}{male_gq} );

            my $marker = $scf_contigs{$scf}{$contig}{markers}{$pos}{marker};
            my $cm =
              defined $map_markers{$marker}{cm}
              ? $map_markers{$marker}{cm}
              : '-';
            if ( $cm ne '-' ) {
                $contig_cm{$cm}++;
            }
            my $cm_order =
              defined $map_markers{$marker}{cm_order}
              ? $map_markers{$marker}{cm_order}
              : '-';

            if ( $cm_order ne '-' ) {
                $contig_cm_order{$cm_order}++;

                if (
                    ( $prev_cm_order ne '-' )
                    and (  ( ( $prev_cm_order + 1 ) < $cm_order )
                        or ( ( $prev_cm_order - 1 ) > $cm_order ) )
                  )
                {
                    if ( $prev_cm_chrpr == $map_markers{$marker}{chrpr} ) {
                        $contig_cm_leaps{$contig}{$cm_order}++;
                        print "BREAK\t$scf\t$cm_order\t$prev_cm_order\n"
                          if $debug;
                    }
                }
                $prev_cm_order = $cm_order;
                $prev_cm_chrpr = $map_markers{$marker}{chrpr};
            }

            my $marker_line =
"$scf\t$contig\t$pos\t$marker\t$cm\t$cm_order\t$scf_contigs{$scf}{$contig}{markers}{$pos}{type}\t$map_markers{$marker}{pos}\t$scf_contigs{$scf}{$contig}{markers}{$pos}{depth}\t$scf_contigs{$scf}{$contig}{markers}{$pos}{qual}\t$scf_contigs{$scf}{$contig}{markers}{$pos}{pos_dp}\t$scf_contigs{$scf}{$contig}{markers}{$pos}{mq}\t$scf_contigs{$scf}{$contig}{markers}{$pos}{male_dp}\t$scf_contigs{$scf}{$contig}{markers}{$pos}{female_dp}\t$sex_dp_diff\t$scf_contigs{$scf}{$contig}{markers}{$pos}{male_gq}\t$scf_contigs{$scf}{$contig}{markers}{$pos}{female_gq}\t$sex_gq_diff\t$map_markers{$marker}{chrpr}\t $map_markers{$marker}{gt}\n";

            $marker_output .= $marker_line;

            print $marker_line if $debug;

        }

        my @contig_cm_orders = sort { $a <=> $b } keys %contig_cm_order;
        for my $i ( 0 .. $#contig_cm_orders - 1 ) {
            if ( $contig_cm_orders[ $i + 1 ] - $contig_cm_orders[$i] > 1 ) {
                $contig_cm_leaps{$contig}{ $contig_cm_orders[$i] }++;
                $contig_cm_leaps{$contig}{ $contig_cm_orders[ $i + 1 ] }++;
            }
        }

        # If more than one chrpr on this contig, there is a contig misassembly
        if ( ( keys %{ $scf_contigs{$scf}{$contig}{chrpr} } ) > 1 ) {

            # Emit previous contigs as split scaffold
            if (   ( $scf_start_contig < $contig )
                && ( $scf_end_contig < $contig ) )
            {
                if ( $scf_end_contig == 0 ) {
                    $rev_scf =
                      emit_scaffold( $scf, $scf_start_contig, $contig - 2,
                        "-", \%split_scf, \%revised_agp, $rev_scf );
                }
                elsif (
                    $scf_chrpr == (
                        sort {
                            $scf_contigs{$scf}{$contig}{chrpr}
                              {$a} <=> $scf_contigs{$scf}{$contig}{chrpr}{$b}
                          }
                          keys %{ $scf_contigs{$scf}{$contig}{chrpr} }
                    )[0]
                  )
                {
                    $rev_scf =
                      emit_scaffold( $scf, $scf_start_contig, $contig - 2,
                        $scf_chrpr, \%split_scf, \%revised_agp, $rev_scf );
                }
                else {
                    $rev_scf =
                      emit_scaffold( $scf, $scf_start_contig, $scf_end_contig,
                        $scf_chrpr, \%split_scf, \%revised_agp, $rev_scf );
                    if ( $scf_end_contig < ( $contig - 2 ) ) {
                        $rev_scf = emit_scaffold(
                            $scf,
                            $scf_end_contig + 2,
                            $contig - 2,
                            "-", \%split_scf, \%revised_agp, $rev_scf
                        );
                    }
                }
            }

            # Emit this contig on its own
            $rev_scf = emit_scaffold( $scf, $contig, $contig, "-", \%split_scf,
                \%revised_agp, $rev_scf );
            $emitted_this_contig++;

            # Move to next contig, starting to build a new scaffold
            $scf_start_contig    = $contig + 2;
            $scf_end_contig      = 0;
            $candidate_scf_chrpr = (
                sort {
                    $scf_contigs{$scf}{$contig}{chrpr}
                      {$a} <=> $scf_contigs{$scf}{$contig}{chrpr}{$b}
                  }
                  keys %{ $scf_contigs{$scf}{$contig}{chrpr} }
            )[-1];    # Last pos chrpr
        }

        # If exactly one chrpr on this contig, could be OK or could vary
        # from previous contig
        elsif ( ( keys %{ $scf_contigs{$scf}{$contig}{chrpr} } ) == 1 ) {

            # Get chromosome print for this contig and set scaffold chrom
            # if not set already
            my $contig_chrpr =
              ( keys %{ $scf_contigs{$scf}{$contig}{chrpr} } )[0];
            if ( $scf_chrpr == 0 ) { $scf_chrpr = $contig_chrpr; }

            # If contig chrpr matches scf chrpr, include this contig in
            # next emitted scaffold
            if (
                ( keys %{ $contig_cm_leaps{$contig} } > 0 )
                or (
                    !(
                           ( $scf_chrpr == $contig_chrpr )
                        or ( $candidate_scf_chrpr == $contig_chrpr )
                    )
                )
              )
            {    # contig chrpr doesn't match scf chrpr; misassembly
                    # Emit previous contigs as split scaffold
                if ( $scf_start_contig < $contig ) {
                    if ( $scf_end_contig == 0 ) {
                        $rev_scf =
                          emit_scaffold( $scf, $scf_start_contig, $contig - 2,
                            "-", \%split_scf, \%revised_agp, $rev_scf );
                    }
                    else {

                        $scf_end_contig =
                          (      ( keys %{ $contig_cm_leaps{$contig} } > 0 )
                              && ( $scf_start_contig < $scf_cm_end_contig ) )
                          ? $scf_cm_end_contig
                          : $scf_end_contig;
                        $rev_scf = emit_scaffold(
                            $scf,
                            $scf_start_contig,
                            $scf_end_contig,
                            (
                                ( $scf_chrpr > 0 )
                                  && (
                                    keys
                                    %{ $contig_cm_leaps{$scf_end_contig} } <=
                                    1 )
                              ) ? $scf_chrpr : '-',
                            \%split_scf,
                            \%revised_agp,
                            $rev_scf
                        );
                        if ( $scf_end_contig < ( $contig - 2 ) ) {
                            $rev_scf = emit_scaffold(
                                $scf,
                                $scf_end_contig + 2,
                                $contig - 2,
                                ( keys %{ $contig_cm_leaps{$contig} } > 0 )
                                ? $scf_chrpr
                                : "-",
                                \%split_scf,
                                \%revised_agp,
                                $rev_scf
                            );
                        }
                    }
                }

                # Start new scaffold
                $scf_start_contig = $contig;
                $scf_end_contig   = $contig;

                # Update scaffold contig
                $scf_chrpr = $contig_chrpr;
            }
            else {
                $scf_end_contig = $contig;
                if ( keys %contig_cm > 0 ) {
                    $scf_cm_end_contig = $contig;
                }
                $scf_chrpr           = $contig_chrpr;
                $candidate_scf_chrpr = 0;
            }
        }
        else {    # No markers for this contig, so no chrpr information
        }

        # If last contig, emit remainder of scaffold
        if (   ( !$emitted_this_contig )
            && ( $contig >= keys %{ $scf_contigs{$scf} } ) )
        {
            if ( $scf_start_contig eq $contig ) {
                $rev_scf = emit_scaffold(
                    $scf,
                    $scf_start_contig,
                    $contig,
                    (
                        (
                            ( keys %{ $scf_contigs{$scf}{$contig}{chrpr} } ) ==
                              1
                        )
                          && ( $scf_chrpr > 0 )
                          && ( keys %{ $contig_cm_leaps{$contig} } <= 1 )
                      )
                    ? $scf_chrpr
                    : '-',
                    \%split_scf,
                    \%revised_agp,
                    $rev_scf
                );

            }
            else {
                $rev_scf = emit_scaffold(
                    $scf,
                    $scf_start_contig,
                    $contig,
                    (
                        ( $scf_chrpr > 0 ) && ( ( $candidate_scf_chrpr == 0 )
                            or ( $candidate_scf_chrpr eq $scf_chrpr ) )
                      )
                    ? $scf_chrpr
                    : '-',
                    \%split_scf,
                    \%revised_agp,
                    $rev_scf
                );
            }
        }
    }

    print "\n" if $debug;
    my $discarded_bases = $scf_lengths{gap}{$scf};
    my $mapped_bases    = 0;
    my $unmapped_bases  = 0;
    if ( ( !defined $split_scf{$scf} ) or ( $split_scf{$scf}{chrpr} eq '-' ) ) {

        printf "$scf (%2d contigs, %7d bp)\n",
          scalar keys %{ $scf_contigs{$scf} }, $scf_lengths{gap}{$scf}
          if $verbose;

        my %split_genes;
        foreach my $scf_name ( sort keys %split_scf ) {
            $discarded_bases -= $split_scf{$scf_name}{length};
            if ( $split_scf{$scf_name}{chrpr} eq '-' ) {
                $unmapped_bases += $split_scf{$scf_name}{length};
            }
            else {
                $mapped_bases += $split_scf{$scf_name}{length};
            }

            foreach my $gene ( keys %{ $scf_genes{$scf} } ) {
                if (
                    (
                        (
                            $scf_genes{$scf}{$gene}{start} <
                            $split_scf{$scf_name}{scf_start}
                        )
                        and ( $scf_genes{$scf}{$gene}{end} >
                            $split_scf{$scf_name}{scf_start} )
                    )
                    or (
                        (
                            $scf_genes{$scf}{$gene}{start} <
                            $split_scf{$scf_name}{scf_end}
                        )
                        && ( $scf_genes{$scf}{$gene}{end} >
                            $split_scf{$scf_name}{scf_end} )
                    )
                  )
                {
                    $split_genes{$gene}++;
                }
            }

            printf
              "\t%-22s (%7d  bp, %7d \- %7d, contigs %2d \- %2d, chr %4s)\n",
              $scf_name,
              $split_scf{$scf_name}{length},
              $split_scf{$scf_name}{scf_start},
              $split_scf{$scf_name}{scf_end},
              $split_scf{$scf_name}{start_ctg},
              $split_scf{$scf_name}{end_ctg}, $split_scf{$scf_name}{chrpr}
              if ($verbose);
        }
        if ($verbose) {
            printf "\t%7d bp mapped\n",   $mapped_bases;
            printf "\t%7d bp unmapped\n", $unmapped_bases;
            if (   ( $discarded_bases < $scf_lengths{gap}{$scf} )
                && ( $discarded_bases > 0 ) )
            {
                printf "\t%7d bp discarded\n", $discarded_bases;
            }
            if ( keys %split_genes > 0 ) {
                print "Broken genes: ";
                foreach my $gene ( sort keys %split_genes ) {
                    print "$gene ";
                }
                print "\n";
            }
            print "\n";
        }
    }
}

map {
    if ( !defined $map_scf_marker{$_} )
    {
        $scf_class{empty}{$_}++;
    }
  }
  keys %all_scf_marker;

map {
    if (   ( !defined $scf_class{rev_unmapped}{$_} )
        && ( !defined $scf_class{rev_mapped}{$_} ) )
    {
        $scf_class{rev_unmapped}{$_}++;
    }
} keys %{ $scf_class{revised} };

foreach my $scf ( keys %scf_contigs ) {
    if (   ( !defined $revised_agp{$scf} )
        && ( !defined $scf_class{mixed}{$scf} ) )
    {
        foreach my $contig ( keys %{ $scf_contigs{$scf} } ) {
            $revised_agp{$scf}{$contig}{scf_start} =
              $scf_contigs{$scf}{$contig}{scf_start};
            $revised_agp{$scf}{$contig}{scf_end} =
              $scf_contigs{$scf}{$contig}{scf_end};
            $revised_agp{$scf}{$contig}{type} =
              $scf_contigs{$scf}{$contig}{type};

            if ( $revised_agp{$scf}{$contig}{type} eq "W" ) {
                $revised_agp{$scf}{$contig}{name} =
                  $scf_contigs{$scf}{$contig}{name};
                $revised_agp{$scf}{$contig}{ctg_start} =
                  $scf_contigs{$scf}{$contig}{ctg_start};
                $revised_agp{$scf}{$contig}{ctg_end} =
                  $scf_contigs{$scf}{$contig}{ctg_end};
                $revised_agp{$scf}{$contig}{orientation} =
                  $scf_contigs{$scf}{$contig}{orientation};
            }
            elsif ( $revised_agp{$scf}{$contig}{type} eq "N" ) {
                $revised_agp{$scf}{$contig}{length} =
                  $scf_contigs{$scf}{$contig}{length};
            }
            $revised_agp{$scf}{$contig}{comment} =
"# $scf\:$scf_contigs{$scf}{$contig}{scf_start}:$scf_contigs{$scf}{$contig}{scf_end}:$contig";
        }
    }
}

sub emit_scaffold {
    my ( $scf, $start, $end, $chrpr, $split_ref, $agp_ref, $rev_scf ) = @_;
    my $scf_name;
    if ( $chrpr eq 846 ) { $chrpr = '-' }
    if (   ( $start == 1 )
        && ( $end == keys %{ $scf_contigs{$scf} } )
        && ( $chrpr ne '-' ) )
    {
        $scf_class{pure}{$scf}++;
        $scf_name = $scf;
        $revised_genome{chrpr}{$scf_name} = $chrpr;
        $scf_class{rev_mapped}{$scf_name}++;
    }
    else {
        if ( ( keys %{ $scf_contigs{$scf} } ) > 1 ) {
            $scf_name = "scf$rev_scf";
            $rev_scf++;
        }
        else {
            $scf_name = $scf;
        }
        if ( defined $scf_class{revised}{$scf} ) {
            delete $scf_class{revised}{$scf};
            delete $revised_genome{length}{$scf};
            delete $revised_genome{length_no_gap}{$scf};
            delete $revised_genome{chrpr}{$scf};
        }
        $scf_class{mixed}{$scf}++;
        $scf_class{added}{$scf_name}++;
        $scf_class{revised}{$scf_name}++;
        if ( $chrpr eq '-' ) {
            $scf_class{added_unmapped}{$scf_name}++;
            $scf_class{rev_unmapped}{$scf_name}++;
        }
        else {
            $scf_class{added_mapped}{$scf_name}++;
            $scf_class{rev_mapped}{$scf_name}++;
        }
    }

    my $scf_length =
      $scf_contigs{$scf}{$end}{scf_end} -
      $scf_contigs{$scf}{$start}{scf_start} + 1;

    # Calculate length without gaps
    my $scf_length_no_gap = 0;
    foreach my $contig ( $start .. $end ) {
        next if ( $scf_contigs{$scf}{$contig}{type} eq "N" );
        $scf_length_no_gap +=
          $scf_contigs{$scf}{$contig}{scf_end} -
          $scf_contigs{$scf}{$contig}{scf_start} + 1;
    }

    $revised_genome{length}{$scf_name}        = $scf_length;
    $revised_genome{length_no_gap}{$scf_name} = $scf_length_no_gap;
    $revised_genome{chrpr}{$scf_name}         = $chrpr;
    $split_ref->{$scf_name}{length}           = $scf_length;
    $split_ref->{$scf_name}{chrpr}            = $chrpr;
    $split_ref->{$scf_name}{start_ctg}        = $start;
    $split_ref->{$scf_name}{end_ctg}          = $end;
    $split_ref->{$scf_name}{scf_start} = $scf_contigs{$scf}{$start}{scf_start};
    $split_ref->{$scf_name}{scf_end}   = $scf_contigs{$scf}{$end}{scf_end};

    my $rev_ctg_num = 1;
    foreach my $contig ( $start .. $end ) {
        $revised_agp{$scf_name}{$rev_ctg_num}{scf_start} =
          $scf_contigs{$scf}{$contig}{scf_start} -
          $scf_contigs{$scf}{$start}{scf_start} + 1;
        $revised_agp{$scf_name}{$rev_ctg_num}{scf_end} =
          $scf_contigs{$scf}{$contig}{scf_end} -
          $scf_contigs{$scf}{$start}{scf_start} + 1;
        $revised_agp{$scf_name}{$rev_ctg_num}{type} =
          $scf_contigs{$scf}{$contig}{type};

        if ( $revised_agp{$scf_name}{$rev_ctg_num}{type} eq "W" ) {
            $revised_agp{$scf_name}{$rev_ctg_num}{name} =
              $scf_contigs{$scf}{$contig}{name};
            $revised_agp{$scf_name}{$rev_ctg_num}{ctg_start} =
              $scf_contigs{$scf}{$contig}{ctg_start};
            $revised_agp{$scf_name}{$rev_ctg_num}{ctg_end} =
              $scf_contigs{$scf}{$contig}{ctg_end};
            $revised_agp{$scf_name}{$rev_ctg_num}{orientation} =
              $scf_contigs{$scf}{$contig}{orientation};

            foreach my $pos ( sort { $a <=> $b }
                keys %{ $scf_contigs{$scf}{$contig}{markers} } )
            {
                my $marker = $scf_contigs{$scf}{$contig}{markers}{$pos}{marker};

                if ( $scf_name ne $scf ) {
                    my $new_marker_pos =
                      $pos - $scf_contigs{$scf}{$start}{scf_start} + 1;
                    if ($debug) {
			print "$scf_name\t$rev_ctg_num\t$new_marker_pos\t$marker\t";
                        print defined( $map_markers{$marker}{cm} )
                          ? "$map_markers{$marker}{cm}"
                          : '-';
                        print "\n";
		    }
	        }
                delete $map_markers{$marker}{scf}{$scf}{$contig};
                if ( $revised_genome{chrpr}{$scf_name} ne '-' ) {
                   $map_markers{$marker}{scf}{$scf_name}{$rev_ctg_num}++;
                }
            }
        }
        elsif ( $revised_agp{$scf_name}{$rev_ctg_num}{type} eq "N" ) {
            $revised_agp{$scf_name}{$rev_ctg_num}{length} =
              $scf_contigs{$scf}{$contig}{length};
        }
        $revised_agp{$scf_name}{$rev_ctg_num}{comment} =
"# $scf\:$scf_contigs{$scf}{$contig}{scf_start}:$scf_contigs{$scf}{$contig}{scf_end}:$contig";
        $rev_ctg_num++;
    }
    return $rev_scf;
}

sub stat_line {
    my ( $scf_list, $filter, $scf_len_ref ) = @_;
    my @scf_lengths;
    my $total_len = 0;
    my $n50_len   = 0;
    my $n50;
    map {
        push @scf_lengths, $scf_len_ref->{$_};
        $total_len += $scf_len_ref->{$_};
      }
      keys %{$scf_list};
    foreach my $scf_len ( sort { $b <=> $a } @scf_lengths ) {
        $n50_len += $scf_len;
        if ( $n50_len >= ( $total_len / 2 ) ) { $n50 = $scf_len; last; }
    }
    my $out_line = sprintf "%4d | %9d | %6d | $filter\n",
      scalar @scf_lengths,
      $total_len, $n50;
    return $out_line;
}

print "Genome statistics with gaps between contigs included:\n";
print_genome_stats( \%scf_class, \%scf_lengths, \%revised_genome, 1 );

print "\n\nGenome statistics with gaps between contigs excluded:\n";
print_genome_stats( \%scf_class, \%scf_lengths, \%revised_genome, 0 );

sub print_genome_stats {
    my ( $scf_class_ref, $scf_lengths_ref, $revised_genome_ref, $gapped ) = @_;
    my $scf_lengths;
    my $revised_scf_lengths;
    if ($gapped) {
        $scf_lengths         = $scf_lengths_ref->{gap};
        $revised_scf_lengths = $revised_genome_ref->{length};
    }
    else {
        $scf_lengths         = $scf_lengths_ref->{no_gap};
        $revised_scf_lengths = $revised_genome_ref->{length_no_gap};
    }

    print "Scf  | Bases     | N50    | Filter\n";
    print stat_line( $scf_class{all}, "whole genome", \%{$scf_lengths} );
    print stat_line(
        $scf_class_ref->{nomarkers},
        "scaffolds without candidate markers",
        \%{$scf_lengths}
    );
    print stat_line(
        $scf_class_ref->{markers},
        "scaffolds with candidate marker(s)",
        \%{$scf_lengths}
    );
    print stat_line(
        $scf_class_ref->{empty},
        "scaffolds with unmapped markers",
        \%{$scf_lengths}
    );
    print stat_line( $scf_class_ref->{map}, "scaffolds with   mapped markers",
        \%{$scf_lengths} );
    print stat_line(
        $scf_class_ref->{pure},
        "scaffolds with  1 chromosome print",
        \%{$scf_lengths}
    );
    print stat_line(
        $scf_class_ref->{mixed},
        "scaffolds with >1 chromosome print",
        \%{$scf_lengths}
    );
    print "\n";
    print stat_line(
        $scf_class_ref->{revised},
        "revised genome",
        \%{$revised_scf_lengths}
    );
    print stat_line(
        $scf_class_ref->{added},
        "scaffolds created from scaffolds with >1 chromosome print",
        \%{$revised_scf_lengths}
    );
    print stat_line(
        $scf_class_ref->{added_unmapped},
"scaffolds created from scaffolds with >1 chromosome print and removed from map",
        \%{$revised_scf_lengths}
    );
    print stat_line(
        $scf_class_ref->{added_mapped},
"scaffolds created from scaffolds with >1 chromosome print remaining on map",
        \%{$revised_scf_lengths}
    );
    print stat_line(
        $scf_class_ref->{rev_unmapped},
        "revised scaffolds without candidate markers",
        \%{$revised_scf_lengths}
    );
    print stat_line(
        $scf_class_ref->{rev_mapped},
        "revised scaffolds mapped",
        \%{$revised_scf_lengths}
    );
    print "\n";
}
open my $scf_out, '>', $scaffold_agp_output
  or croak "Can't open scaffold AGP file $scaffold_agp_output! $OS_ERROR";

foreach my $scf ( sort keys %revised_agp ) {
    foreach my $contig ( sort { $a <=> $b } keys %{ $revised_agp{$scf} } ) {
        print $scf_out
"$scf\t$revised_agp{$scf}{$contig}{scf_start}\t$revised_agp{$scf}{$contig}{scf_end}\t$contig\t";
        if ( $revised_agp{$scf}{$contig}{type} eq "W" ) {
            print $scf_out
"W\t$revised_agp{$scf}{$contig}{name}\t$revised_agp{$scf}{$contig}{ctg_start}\t$revised_agp{$scf}{$contig}{ctg_end}\t$revised_agp{$scf}{$contig}{orientation}\t$revised_agp{$scf}{$contig}{comment}\n";
        }
        elsif ( $revised_agp{$scf}{$contig}{type} eq "N" ) {
            print $scf_out
"N\t$revised_agp{$scf}{$contig}{length}\tfragment\tyes\t\t$revised_agp{$scf}{$contig}{comment}\n";
        }
    }
}
close $scf_out;

open my $lg_file, '<', $lg_filename
  or croak "Can't open $lg_filename $OS_ERROR!\n";
my %chrom_scf;
while ( my $chrom_line = <$lg_file> ) {
    chomp $chrom_line;
    my @fields = split /,/, $chrom_line;
    if (   ( defined $chrom_scf{ $fields[2] }{chrom} )
        && ( $chrom_scf{ $fields[2] }{chrom} ne $fields[0] ) )
    {
        print STDERR
"Scaffold $fields[2] linked to more than one linkage group: $fields[0] $chrom_scf{$fields[2]}{chrom}\n";
        exit;
    }
    $chrom_scf{ $fields[2] }{chrom} = $fields[0];
}
close $lg_file;

my %ll_chrom;
my %chrom_ll;
foreach my $scf ( keys %map_scf_marker ) {
    next if ( defined $scf_class{mixed}{$scf} );
    if ( defined $chrom_scf{$scf} ) {
        my $chrom = $chrom_scf{$scf}{chrom};
        my $chrpr;
        foreach my $marker ( keys %{ $map_scf_marker{$scf} } ) {
            if (   ( defined $chrpr )
                && ( $chrpr ne $map_markers{$marker}{chrpr} ) )
            {
                print STDERR
"Scaffold $scf has conflicting chromosome prints $chrpr $map_markers{$marker}{chrpr}\n";
            }
            $chrpr = $map_markers{$marker}{chrpr};
        }
        $ll_chrom{$chrpr}{$chrom}{$scf}++;
        $chrom_ll{$chrom}{$chrpr}{$scf}++;
    }
}

foreach my $ll_marker ( sort { $a <=> $b } keys %ll_chrom ) {
    print "$ll_marker" if $debug;
    foreach my $chrom (
        reverse
        sort { $ll_chrom{$ll_marker}{$a} <=> $ll_chrom{$ll_marker}{$b} }
        keys %{ $ll_chrom{$ll_marker} }
      )
    {
        print "\t$chrom:" . ( scalar keys %{ $ll_chrom{$ll_marker}{$chrom} } )
          if $debug;
    }
    print "\n" if $debug;
}
print "\n\n" if $debug;
foreach my $chrom ( sort { $a <=> $b } keys %chrom_ll ) {
    foreach my $ll_marker (
        reverse sort { $chrom_ll{$chrom}{$a} <=> $chrom_ll{$chrom}{$b} }
        keys %{ $chrom_ll{$chrom} }
      )
    {
        foreach my $scf ( keys %{ $chrom_ll{$chrom}{$ll_marker} } ) {
            print "$chrom\t$ll_marker\t$scf\n" if $debug;
        }
    }
    print "\n" if $debug;
}

# Construct chromosome AGP

open my $chrom_out, '>', $chromosome_agp_output
  or croak "Can't open chromosome AGP file $chromosome_agp_output! $OS_ERROR\n";

print "   \t  \tGROUPED ON CHROMOSOME   \t        \tLINKED ON CHROMOSOME\n";
print
"Chr\tcM\tScaffolds\tGapped Bases\tUngapped Bases\tScaffolds\tGapped Bases\tUngapped Bases\n";
my $total_scf                     = 0;
my $total_bases                   = 0;
my $total_bases_ungapped          = 0;
my $total_cm                      = 0;
my $total_unmapped_scf            = 0;
my $total_unmapped_bases          = 0;
my $total_unmapped_bases_ungapped = 0;
foreach my $chr ( sort { $a <=> $b } keys %chrom_ll ) {
    my $chrpr = ( keys %{ $chrom_ll{$chr} } )[0];
    my $lg    = $chrpr_to_lg{$chrpr};
    my %chrom_agp;
    my %scf_ctg_order;
    my @cm_list = sort { $a <=> $b } keys %{ $map{$lg} };
    my @chr_scf_order;
    my %scf_comments;
    my %scf_comment_contigs;
    my $last_cm         = 0;
    my $chrom_scf_count = 0;

    foreach my $cm (@cm_list) {
        foreach my $id ( keys %{ $map{$lg}{$cm} } ) {
            foreach my $scf ( keys %{ $map_markers{$id}{scf} } ) {

                foreach my $contig (
                    sort { $a <=> $b }
                    keys %{ $map_markers{$id}{scf}{$scf} }
                  )
                {
                    $chrom_agp{$cm}{scf}{$scf}{$contig}++;
                    $scf_ctg_order{$scf}{$contig}{$cm} = $map_markers{$id}{pos};
                    $scf_comment_contigs{$scf}{$cm}{$contig}++;
                    $last_cm = $cm;
                }
            }
        }
    }

    foreach my $scf ( keys %scf_comment_contigs ) {
        $scf_comments{$scf} = " Contigs |";
        foreach
          my $cm ( sort { $a <=> $b } keys %{ $scf_comment_contigs{$scf} } )
        {
            $scf_comments{$scf} .= " $cm cM:";
            foreach my $contig (
                sort { $a <=> $b }
                keys %{ $scf_comment_contigs{$scf}{$cm} }
              )
            {
                $scf_comments{$scf} .= " $contig";
            }
            $scf_comments{$scf} .= " |";
        }
    }

    my %scf_orient;

    foreach my $scf ( keys %scf_ctg_order ) {
        my $cm_start;
      CONTIG:
        foreach my $contig ( sort { $a <=> $b } keys %{ $scf_ctg_order{$scf} } )
        {
            foreach my $cm (
                sort { $scf_ctg_order{$scf}{$contig}{$a} <=> $scf_ctg_order{$scf}{$contig}{$b} }
                keys %{ $scf_ctg_order{$scf}{$contig} }
              )
            {
                if ( !defined $cm_start ) {
                    $cm_start = $cm;
                    next CONTIG;
                }
                if ( $cm_start < $cm ) {
                    $scf_orient{$scf} = '+';
                }
                if ( $cm_start > $cm ) {
                    $scf_orient{$scf} = '-';
                }
            }

        }
        if ( !defined $scf_orient{$scf} ) {
            $scf_orient{$scf} = "+";
        }
    }

    my %cm_scf_order;
    for my $i ( 0 .. ($#cm_list) ) {
        my @last_scf;
        foreach my $scf ( keys %{ $chrom_agp{ $cm_list[$i] }{scf} } ) {
            if (   ( $i ne $#cm_list )
                && ( defined $chrom_agp{ $cm_list[ $i + 1 ] }{scf}{$scf} ) )
            {

                $chrom_agp{ $cm_list[ $i + 1 ] }{ignore}{$scf}++;

                if ( !defined $chrom_agp{ $cm_list[$i] }{ignore}{$scf} ) {
                    push @last_scf, $scf;
                }
            }
            else {
                if ( !defined $chrom_agp{ $cm_list[$i] }{ignore}{$scf} ) {
                    push @{ $cm_scf_order{ $cm_list[$i] } }, $scf;
                    push @chr_scf_order, $scf;
                }
            }
        }
        if ( @last_scf > 0 ) {
            push @{ $cm_scf_order{ $cm_list[$i] } }, @last_scf;
            push @chr_scf_order, @last_scf;
        }
    }

    my $chr_loc        = 1;
    my $component_num  = 0;
    my $prev_cm        = 0;
    my $bases_gapped   = 0;
    my $bases_ungapped = 0;

    foreach my $cm ( sort { $a <=> $b } keys %chrom_agp ) {
        my $num_cm_scfs =
          defined $cm_scf_order{$cm} ? @{ $cm_scf_order{$cm} } : 0;
        if (   ( $component_num > 0 )
            && ( @chr_scf_order > 0 )
            && ( $num_cm_scfs > 0 ) )
        {

            # Add gap between cM appropriate to distance between
            # current cM and previous cM
            print $chrom_out "chr$chr\t$chr_loc";
            $chr_loc += 99;
            print $chrom_out "\t$chr_loc";
            $chr_loc++;
            $component_num++;
            print $chrom_out "\t$component_num\tN\t100\tfragment\tno";
            printf $chrom_out
"\t# %6.3f to %6.3f cM ------------------------------------------------------",
              $prev_cm, $cm;
            print $chrom_out "\n";
        }

        my $scf_num = 0;
        my $last_scf;

        foreach my $scf ( @{ $cm_scf_order{$cm} } ) {
            my $agp_scf;
            my $agp_gap;

            # Output scaffold
            $agp_scf .= "chr$chr\t$chr_loc";
            $chr_loc += $revised_genome{length}{$scf} - 1;
            $agp_scf .= "\t$chr_loc";
            $chr_loc++;
            $component_num++;
            $agp_scf .=
              "\t$component_num\tD\t$scf\t1\t$revised_genome{length}{$scf}";
            $agp_scf .= "\t$scf_orient{$scf}\t#$scf_comments{$scf}";
            $agp_scf .= "\n";
            $scf_num++;
            $chrom_scf_count++;
            $bases_ungapped += $revised_genome{length_no_gap}{$scf};
            print $chrom_out $agp_scf;
            shift @chr_scf_order;

            # Output gap if not the last scaffold at this cM
            if ( $scf_num < $num_cm_scfs ) {
                $agp_gap = "chr$chr\t$chr_loc";
                $chr_loc += 99;
                $agp_gap .= "\t$chr_loc";
                $chr_loc++;
                $component_num++;
                $agp_gap .= "\t$component_num\tN\t100\tfragment\tno\t# $cm cM";
                $agp_gap .= "\n";
                print $chrom_out $agp_gap;
            }
        }
        $prev_cm = $cm;
    }

    $bases_gapped = $chr_loc - 1;
    printf "%3d\t%7.1f\t%9d\t%9d\t%9d", $chr, $last_cm, $chrom_scf_count,
      $bases_gapped, $bases_ungapped;

    $total_scf            += $chrom_scf_count;
    $total_bases          += $bases_gapped;
    $total_cm             += $last_cm;
    $total_bases_ungapped += $bases_ungapped;

    $component_num = 0;
    $chr_loc       = 1;
    my $gap                     = "";
    my $unmapped_scf            = 0;
    my $unmapped_bases          = 0;
    my $unmapped_bases_ungapped = 0;
    foreach my $scf ( keys %{ $revised_genome{chrpr} } ) {
        if (   ( $chrpr eq $revised_genome{chrpr}{$scf} )
            && ( !defined $scf_ctg_order{$scf} ) )
        {
            if ( $gap ne "" ) { print $chrom_out $gap; $gap = ""; }
            print $chrom_out "chr$chr\_unmapped\t$chr_loc";
            $chr_loc += $revised_genome{length}{$scf} - 1;
            print $chrom_out "\t$chr_loc";
            $chr_loc++;
            $component_num++;
            print $chrom_out
              "\t$component_num\tD\t$scf\t1\t$revised_genome{length}{$scf}";
            print $chrom_out "\t+";
            print $chrom_out "\n";

            $unmapped_scf++;
            $unmapped_bases          += $revised_genome{length}{$scf};
            $unmapped_bases_ungapped += $revised_genome{length_no_gap}{$scf};
            $gap = "chr$chr\_unmapped\t$chr_loc";
            $chr_loc += 99;
            $gap .= "\t$chr_loc";
            $chr_loc++;
            $component_num++;
            $gap .= "\t$component_num\tN\t100\tfragment\tno\t";
            $gap .= "\n";
        }
    }

    printf "\t%9d\t%9d\t%9d\n", $unmapped_scf, $unmapped_bases,
      $unmapped_bases_ungapped;
    $total_unmapped_scf            += $unmapped_scf;
    $total_unmapped_bases          += $unmapped_bases;
    $total_unmapped_bases_ungapped += $unmapped_bases_ungapped;
}

printf "ALL\t%7.1f\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\n", $total_cm, $total_scf,
  $total_bases, $total_bases_ungapped,
  $total_unmapped_scf, $total_unmapped_bases, $total_unmapped_bases_ungapped;

print "Total scaffolds        = " . ( $total_scf + $total_unmapped_scf ) . "\n";
print "Total gapped bases     = "
  . ( $total_bases + $total_unmapped_bases ) . "\n";
print "Total ungapped bases   = "
  . ( $total_bases_ungapped + $total_unmapped_bases_ungapped ) . "\n";
