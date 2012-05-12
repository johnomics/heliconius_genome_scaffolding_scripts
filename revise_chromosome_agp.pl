#!/usr/bin/env perl

# revise_chromosome_agp.pl
#
# This script is provided FOR INFORMATION ONLY and is not intended for use.
# It was written specifically for the Heliconius melpomene genome paper
# (Heliconius Genome Consortium, doi: 10.1038/nature11041) and has not
# been adapted for general use.

# Purpose: Take whole genome chromosome AGP file
#          and individual Bambus chromosome AGPs
#          and output revised AGP file based on new Bambus links

# Input  : chromosome AGP file
#          directory containing Bambus2 AGP files
#          Z chromosome AGP file

# Author: John Davey john.davey@ed.ac.uk
# Begun 5/12/11

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
use List::Util qw/min max/;

# Autoflush output so reporting on progress works
$| = 1;

my $chrom_agp_filename = "";
my $bambus_dir         = "";
my $z_agp_filename     = "";

my $options_okay = GetOptions(
    'chrom_agp=s'  => \$chrom_agp_filename,
    'bambus_dir=s' => \$bambus_dir,
    'z_agp=s'      => \$z_agp_filename,
);

croak
  "\nUsage: perl revise_chromosome_agp.pl -c chrom_agp_filename -b bambus_dir\n"
  if !$options_okay;

croak "No chromosome AGP file! Please specify -c $OS_ERROR\n"
  if ( $chrom_agp_filename eq '' );
croak "No Bambus output directory! Please specify -b $OS_ERROR\n"
  if ( $bambus_dir eq '' );
croak "No Z chromosome AGP file! Please specify -z $OS_ERROR\n" if ($z_agp_filename eq '');


print "# CHROMOSOME AGP GENERATED ON " . localtime() . "\n";
print "#\n";
print "# This is not a conventional AGP file. It describes an ordering of scaffolds\n";
print "# based on RAD linkage map data. Not all scaffolds are ordered and attention\n";
print "# should be paid to these comments and the comments associated with each\n";
print "# scaffold and gap.\n#\n";
print "# Each scaffold is assigned to a chromosome. Most scaffolds are located on\n";
print "# the RAD linkage map and can be assigned to one or more cM loci, shown in\n";
print "# the scaffold comment. The comment shows the contigs of the scaffold that\n";
print "# map to each cM locus. Groups of scaffolds at different cM loci can be\n";
print "# considered to be reliably ordered, but multiple scaffolds at one cM locus\n";
print "# may not be. An attempt has been made to order these scaffolds using mate\n";
print "# pair data.\n#\n";
print "# Each gap comment describes the connection between the neighbouring scaffolds.\n";
print "# The possible gap comments, in order of strength from weak to strong, are:\n";
print "# recombination (400bp) - the two scaffolds are separated by a recombination.\n";
print "#                         They may be physically next to each other but there\n";
print "#                         is no evidence to decide this.\n";
print "# unordered     (300bp) - the two scaffolds are at the same cM locus, but\n";
print "#                         there is no evidence to order them within this locus.\n";
print "# ordered       (200bp) - there is only one possible ordering of these scaffolds,\n";
print "#                         for example if there is only one scaffold at a\n";
print "#                         particular cM locus. But there is no direct connecting\n";
print "#                         evidence and it may be that other scaffolds should be\n";
print "#                         inserted between these two scaffolds.\n";
print "# bridged       (100bp) - the two scaffolds are connected by mate pair data and\n";
print "#                         so can be confidently bridged together.\n#\n";
print "# Some chromosomes have an additional set of 'unmapped' scaffolds, listed under\n";
print "# chrN_unmapped. These scaffolds are linked to the chromosome by a chromosome\n";
print "# print marker but have no other markers positioning them on the chromosome.\n";
print "# No ordering can be inferred or assumed for these scaffolds.\n";
print "#\n";


my %gap_size = (
    'bridged'       => 100,
    'ordered'       => 200,
    'unordered'     => 300,
    'recombination' => 400
);

open my $chrom_agp, '<', $chrom_agp_filename
  or croak "Can't open $chrom_agp_filename: $OS_ERROR!\n";

my %chr_cMs;
my %scf_cMs;
my %chr_cMs_scf;
my %scf_orient;
my %scf_comment;
my %scf_lengths;
my %chr_cM_scfpos;
my %unmapped_scfs;

while ( my $chr_agp_line = <$chrom_agp> ) {
    chomp $chr_agp_line;
    my @agp_fields = split /\t/, $chr_agp_line;
    next if ( $agp_fields[4] eq 'N' );
    my $chr = $agp_fields[0];
    my $scf = $agp_fields[5];
    $scf_orient{$scf}  = $agp_fields[8];
    $scf_comment{$scf} = $agp_fields[9];
    $scf_lengths{$scf} = $agp_fields[7];

    my @cMs;
    if ( $chr !~ /unmapped/ ) {
        my @loci = split /\|/, $agp_fields[9];
        shift @loci;    # Remove '# Contigs :'
        foreach my $locus (@loci) {
            if ( $locus =~ /\ (.+)\ cM:\ (.+)\ / ) {
                my $cur_cM = $1;
                push @cMs, $cur_cM;
                if ( !defined $chr_cM_scfpos{$chr}{$cur_cM} ) {
                    $chr_cM_scfpos{$chr}{$cur_cM} = 1;
                }
            }
        }
    }
    else {
        my $chr_unmapped = $chr;
        $chr_unmapped =~ s/_unmapped//;
        $unmapped_scfs{$chr_unmapped}{$scf}++;
    }

    # If scaffold spans more than one cM, it can be anchored

    foreach my $cM (@cMs) {
        my $pos =
          ( @cMs > 1 )
          ? $chr_cM_scfpos{$chr}{$cM}
          : $chr_cM_scfpos{$chr}{$cM} * -1;
        $chr_cMs{$chr}{$cM}{$pos}     = $scf;
        $chr_cMs_scf{$chr}{$cM}{$scf} = $pos;

        $chr_cM_scfpos{$chr}{$cM}++;
        $scf_cMs{$scf}{$cM}++;
    }
}

close $chrom_agp;

# Load Bambus data

my %scf_bridges;
my %scf_bambus_orient;
foreach my $chr ( sort keys %chr_cMs ) {

    my $scf1      = '';
    my $scf2      = '';
    my $cur_group = 0;
    open my $bambus_agp, '<', "$bambus_dir/$chr.agp"
      or croak "Can't open $bambus_dir/$chr.agp $OS_ERROR\n";
    my %chrom_scf_orient;
    while ( my $bambus_line = <$bambus_agp> ) {
        next if ( $bambus_line =~ /^#/ );
        chomp $bambus_line;
        my @f = split /\t/, $bambus_line;

        my $group = $f[0];
        if ( $cur_group != $group ) {
            $cur_group = $group;
            $scf1      = '';
            $scf2      = '';
        }

        my $scf = $f[5];
        if ( $scf eq 'fragment' ) {
            my $size = $f[4];
            if ( abs($size) > 8000 ) {
                $scf1 = '';
                <$bambus_agp>;
            }
        }
        else {
            $chrom_scf_orient{$scf} = $f[8];
            if ( $scf1 eq '' ) {
                $scf1 = $scf;
            }
            else {
                $scf2 = $scf;
            }
        }

        if ( ( $scf1 ne '' ) && ( $scf2 ne '' ) ) {

            # Record bridge between scf1 and scf2
            $scf_bridges{$scf1}{$scf2}++;
            $scf_bridges{$scf2}{$scf1}++;
            $scf_bambus_orient{$scf1} = $chrom_scf_orient{$scf1};
            $scf_bambus_orient{$scf2} = $chrom_scf_orient{$scf2};
            $scf1                     = $scf2;
            $scf2                     = '';
        }
    }
    close $bambus_agp;
}

my %scf_output;
foreach my $chr (
    sort { ( substr $a, 3 ) <=> ( substr $b, 3 ) }
    keys %chr_cMs
  )
{

    my %bridged_scf;
    my %chr_loc = ( pos => 1, component => 1 );
    foreach my $cM ( sort { $a <=> $b } keys %{ $chr_cMs{$chr} } ) {
        my %start_pos;
        my %unanchored_pos;
        my %end_pos;

        # Split positions into anchored (start and end) or unanchored
        foreach my $pos (
            sort { $a <=> $b }
            keys %{ $chr_cMs{$chr}{$cM} }
          )
        {
            if ( $pos < 0 ) { $unanchored_pos{$pos}++; }
            else {
                my $scf = $chr_cMs{$chr}{$cM}{$pos};
                $end_pos{$pos}++
                  if ( $cM == min( keys %{ $scf_cMs{$scf} } ) );
                $start_pos{$pos}++
                  if ( $cM == max( keys %{ $scf_cMs{$scf} } ) );
            }
        }

        # Walk out from start scaffold for this cM
        if ( keys %start_pos > 0 ) {
            anchor_unanchored_scaffolds(
                \%start_pos,             1,
                $chr_cMs{$chr}{$cM},     $cM,
                $chr_cMs_scf{$chr}{$cM}, \%scf_bridges,
                \%unanchored_pos,        \%bridged_scf,
                \%scf_lengths,           $unmapped_scfs{$chr},
            );
        }

        # Walk in from end scaffold for this cM
        if ( keys %end_pos > 0 ) {
            anchor_unanchored_scaffolds(
                \%end_pos,               -1,
                $chr_cMs{$chr}{$cM},     $cM,
                $chr_cMs_scf{$chr}{$cM}, \%scf_bridges,
                \%unanchored_pos,        \%bridged_scf,
                \%scf_lengths,           $unmapped_scfs{$chr},
            );
        }

 # If no unanchored scaffolds left, try anchoring across start and end scaffolds
        if ( keys %unanchored_pos == 0 ) {
            bridge_anchored_scaffolds( \%start_pos, \%end_pos,
                $chr_cMs{$chr}{$cM},
                \%scf_bridges );
        }

        # Output scaffolds anchored from start of cM
        foreach my $pos ( sort { $a <=> $b } keys %start_pos ) {
            my $scf = $chr_cMs{$chr}{$cM}{$pos};
            my $gap_type =
                ( $chr_loc{component} == 1 )   ? 'ignore'
              : ( defined $bridged_scf{$scf} ) ? 'bridged'
              :                                  'unordered';
            output_scf(
                $scf,              $gap_type,          $chr,
                $cM,               \%chr_loc,          $scf_lengths{$scf},
                $scf_orient{$scf}, $scf_comment{$scf}, \%scf_output,
                \%scf_bridges
            );
        }
        my $original_unanchored_pos_num = keys %unanchored_pos;

        # Output unanchored scaffolds
        bridge_and_output_unanchored_scaffolds(
            \%unanchored_pos,        $chr_cMs{$chr}{$cM},
            $chr_cMs_scf{$chr}{$cM}, $chr,
            $cM, ( scalar keys %start_pos ),
            \%chr_loc,       \%scf_lengths,
            \%scf_orient,    \%scf_comment,
            \%scf_output,    \%scf_bridges,
            \%unmapped_scfs, \%scf_cMs
        );

        # Output scaffolds anchored from end of cM
        my $first_end_pos = 1;
        foreach my $pos ( sort { $a <=> $b } keys %end_pos ) {
            my $scf = $chr_cMs{$chr}{$cM}{$pos};
            my $gap_type =
              ( $chr_loc{component} == 1 ) ? 'ignore'
              : (    ($first_end_pos)
                  && ( keys %start_pos == 0 )
                  && ( $original_unanchored_pos_num == 0 ) ) ? 'recombination'
              : ( ($first_end_pos) && ( $original_unanchored_pos_num == 1 ) )
              ? 'ordered'
              : (    ($first_end_pos)
                  && ( keys %start_pos > 0 )
                  && ( $original_unanchored_pos_num == 0 ) ) ? 'ordered'
              : ( defined $bridged_scf{$scf} ) ? 'bridged'
              :                                  'unordered';
            output_scf(
                $scf,              $gap_type,          $chr,
                $cM,               \%chr_loc,          $scf_lengths{$scf},
                $scf_orient{$scf}, $scf_comment{$scf}, \%scf_output,
                \%scf_bridges
            );
        }
    }

    %chr_loc = ( pos => 1, component => 1 );

    my $output_gap = 0;
    foreach my $unmapped_scf ( keys %{ $unmapped_scfs{$chr} } ) {
        if ($output_gap) {
            print "$chr\_unmapped\t$chr_loc{pos}";
            $chr_loc{pos} += 99;
            print
              "\t$chr_loc{pos}\t$chr_loc{component}\tN\t100\tfragment\tno\t\n";
            $chr_loc{pos}++;
            $chr_loc{component}++;
        }

        print "$chr\_unmapped\t$chr_loc{pos}";
        $chr_loc{pos} += $scf_lengths{$unmapped_scf} - 1;
        print
"\t$chr_loc{pos}\t$chr_loc{component}\tD\t$unmapped_scf\t1\t$scf_lengths{$unmapped_scf}\t+\n";
        $chr_loc{pos}++;
        $chr_loc{component}++;
        $output_gap++;
    }

}

open my $z_agp, '<', $z_agp_filename or croak "Can't open $z_agp_filename: $OS_ERROR\n";
print <$z_agp>;
close $z_agp;


sub bridge_and_output_unanchored_scaffolds {
    my (
        $unanchored_pos_ref, $chr_cM_ref,     $chr_cMs_scf_ref,
        $chr,                $cM,             $start_pos_num,
        $chr_loc_ref,        $scf_length_ref, $scf_orient_ref,
        $scf_comment_ref,    $scf_output_ref, $scf_bridge_ref,
        $unmapped_scfs_ref,  $scf_cMs_ref,
    ) = @_;

    # Check for bridges in unanchored scaffolds

    my $orig_unanchored_pos_num = keys %{$unanchored_pos_ref};

    my %seen_pos;
    my %bridged_scf_links;
    foreach my $pos ( sort { $b <=> $a } keys %{$unanchored_pos_ref} ) {
        my @bridged_scfs;

        my $bridged_scf_num = 0;

        push @bridged_scfs, $chr_cM_ref->{$pos};
        while ( @bridged_scfs > $bridged_scf_num ) {
            $bridged_scf_num = @bridged_scfs;
            foreach my $scf (@bridged_scfs) {
                my $scf_pos = $chr_cMs_scf_ref->{$scf};
                foreach my $scf2 ( keys %{ $scf_bridge_ref->{$scf} } ) {
                    my $scf2_pos = $chr_cMs_scf_ref->{$scf2};
                    next
                      if ( ( defined $scf2_pos )
                        && ( defined $seen_pos{$scf2_pos} ) );

                    $bridged_scf_links{$scf}{$scf2}++;
                    $bridged_scf_links{$scf2}{$scf}++;

                    $scf_orient{$scf}  = $scf_bambus_orient{$scf};
                    $scf_orient{$scf2} = $scf_bambus_orient{$scf2};

                    if ( defined $scf_pos ) {
                        delete $unanchored_pos_ref->{$scf_pos};
                    }
                    if ( defined $scf2_pos ) {
                        delete $unanchored_pos_ref->{$scf2_pos};
                    }
                }
            }
            @bridged_scfs = keys %bridged_scf_links;
            $seen_pos{$pos}++;
        }
    }

    my %seen_scf;
    my $first_unanchored_pos = 1;

    my @last_chain;
    while ( keys %bridged_scf_links > 0 ) {
        my @scf_chain;
        my $scf;
        foreach my $bridged_scf ( keys %bridged_scf_links ) {
            if ( keys %{ $bridged_scf_links{$bridged_scf} } == 1 ) {
                $scf = $bridged_scf;
                last;
            }
        }
        while ( defined $scf ) {
            if ( !$seen_scf{$scf} ) { push @scf_chain, $scf; }
            if ( !defined $bridged_scf_links{$scf} ) { $scf = undef; next; }
            my $scf2;
            while ( !defined $scf2 ) {
                last if ( keys %{ $bridged_scf_links{$scf} } == 0 );
                $scf2 = ( keys %{ $bridged_scf_links{$scf} } )[0];
                if ( $seen_scf{$scf2} ) {
                    $scf2 = undef;
                }
            }
            if ( defined $scf2 ) {
                push @scf_chain, $scf2;
                if ( defined $unmapped_scfs_ref->{$chr}{$scf} ) {
                    delete $unmapped_scfs_ref->{$chr}{$scf};
                }

                if ( defined $unmapped_scfs_ref->{$chr}{$scf2} ) {
                    delete $unmapped_scfs_ref->{$chr}{$scf2};
                }
                delete $bridged_scf_links{$scf}{$scf2};
                delete $bridged_scf_links{$scf2}{$scf};
                if ( keys %{ $bridged_scf_links{$scf} } == 0 ) {
                    delete $bridged_scf_links{$scf};
                }
                if ( keys %{ $bridged_scf_links{$scf2} } == 0 ) {
                    delete $bridged_scf_links{$scf2};
                }
                $seen_scf{$scf}++;
                $seen_scf{$scf2}++;
                $scf = $scf2;
            }
            else {
                $scf = undef;
                if ( keys %{ $bridged_scf_links{$scf} } == 0 ) {
                    delete $bridged_scf_links{$scf};
                }
            }
            if ( keys %{ $bridged_scf_links{$scf} } == 0 ) {
                delete $bridged_scf_links{$scf};
            }
        }

        my %chain_cMs;
        foreach my $scf (@scf_chain) {
            map { $chain_cMs{$_}{$scf}++ } keys %{ $scf_cMs{$scf} };
        }
        if ( keys %chain_cMs > 1 ) {
            if ( $scf_cMs{ $scf_chain[0] } eq max( keys %chain_cMs ) ) {
                @last_chain = reverse @scf_chain;
            }
            else {
                @last_chain = @scf_chain;
            }
            next;
        }

        foreach my $scf (@scf_chain) {
            my $gap_type =
                ( $chr_loc_ref->{component} == 1 ) ? 'ignore'
              : ( $first_unanchored_pos && ( $start_pos_num == 0 ) )
              ? 'recombination'
              : ( $scf_chain[0] eq $scf ) ? 'unordered'
              :                             'bridged';
            $first_unanchored_pos = 0;

            output_scf(
                $scf,                    $gap_type,
                $chr,                    $cM,
                $chr_loc_ref,            $scf_length_ref->{$scf},
                $scf_orient_ref->{$scf}, $scf_comment_ref->{$scf},
                $scf_output_ref,         $scf_bridge_ref
            );
        }
    }

    foreach my $pos ( sort { $b <=> $a } keys %{$unanchored_pos_ref} ) {
        my $scf = $chr_cM_ref->{$pos};
        my $gap_type =
            ( $chr_loc_ref->{component} == 1 ) ? 'ignore'
          : ( $first_unanchored_pos && ( $start_pos_num == 0 ) )
          ? 'recombination'
          : ( $orig_unanchored_pos_num == 1 ) ? 'ordered'
          :                                     'unordered';
        $first_unanchored_pos = 0;
        output_scf(
            $scf,                    $gap_type,
            $chr,                    $cM,
            $chr_loc_ref,            $scf_length_ref->{$scf},
            $scf_orient_ref->{$scf}, $scf_comment_ref->{$scf},
            $scf_output_ref,         $scf_bridge_ref
        );
    }

    foreach my $scf (@last_chain) {
        my $gap_type =
            ( $chr_loc_ref->{component} == 1 ) ? 'ignore'
          : ( $first_unanchored_pos && ( $start_pos_num == 0 ) )
          ? 'recombination'
          : ( $last_chain[0] eq $scf ) ? 'unordered'
          :                              'bridged';
        $first_unanchored_pos = 0;

        output_scf(
            $scf,                    $gap_type,
            $chr,                    $cM,
            $chr_loc_ref,            $scf_length_ref->{$scf},
            $scf_orient_ref->{$scf}, $scf_comment_ref->{$scf},
            $scf_output_ref,         $scf_bridge_ref
        );
    }

    return;
}

sub bridge_anchored_scaffolds {
    my ( $start_pos_ref, $end_pos_ref, $cM_ref, $scf_bridges_ref,
        $bridged_scf_ref )
      = @_;
    my $start_bridge_pos = max( keys %{$start_pos_ref} );
    my $end_bridge_pos   = min( keys %{$end_pos_ref} );
    if (   ( defined $start_bridge_pos )
        && ( defined $end_bridge_pos ) )
    {
        my $start_bridge_scf = $cM_ref->{$start_bridge_pos};
        my $end_bridge_scf   = $cM_ref->{$end_bridge_pos};

        if ( defined $scf_bridges_ref->{$start_bridge_scf}{$end_bridge_scf} ) {
            $bridged_scf_ref->{$end_bridge_scf}++;
        }
    }
    return;
}

sub anchor_unanchored_scaffolds {
    my (
        $walk_pos_ref,   $dir,             $cM_ref,
        $cM,             $scf_cM_ref,      $scf_bridges_ref,
        $unanchored_ref, $bridged_scf_ref, $scf_length_ref,
        $unmapped_scf_ref
    ) = @_;

    my $pos =
      ( $dir == 1 )
      ? min( keys %{$walk_pos_ref} )
      : max( keys %{$walk_pos_ref} );

    my $added_scf = 1;
    while ($added_scf) {
        $added_scf = 0;
        my $anchored_scf = $cM_ref->{$pos};
        if ( defined $scf_bridges_ref->{$anchored_scf} ) {
            foreach my $scf2 ( keys %{ $scf_bridges_ref->{$anchored_scf} } ) {
                if ( defined $scf_cM_ref->{$scf2} ) {
                    my $scf2_pos = $scf_cM_ref->{$scf2};

                    if ( defined( $unanchored_ref->{$scf2_pos} ) ) {
                        my $new_pos = $pos + $dir;
                        $cM_ref->{$new_pos}  = $scf2;
                        $scf_cM_ref->{$scf2} = $new_pos;
                        if ( $dir == 1 ) { $bridged_scf_ref->{$scf2}++ }
                        else             { $bridged_scf_ref->{$anchored_scf}++ }

                        $scf_orient{$scf2} = get_unanchored_orientation(
                            $scf_orient{$anchored_scf},
                            $scf_bambus_orient{$anchored_scf},
                            $scf_bambus_orient{$scf2}
                        );

                        my $scf_to_swap = $cM_ref->{ $new_pos * -1 };
                        delete $cM_ref->{ $new_pos * -1 };

                        if ( $scf_to_swap ne $scf2 ) {
                            $cM_ref->{$scf2_pos}        = $scf_to_swap;
                            $scf_cM_ref->{$scf_to_swap} = $scf2_pos;
                        }
                        delete $unanchored_ref->{ $new_pos * -1 };
                        $walk_pos_ref->{$new_pos}++;

                        $added_scf = 1;
                        $pos       = $new_pos;
                    }
                }
                elsif ( defined $unmapped_scf_ref->{$scf2} ) {
                    # Ignore
                }
            }
        }
    }

    return;
}

sub get_unanchored_orientation {
    my ( $an_o, $an_bo, $unan_bo ) = @_;

    return ( ( $an_o eq '+' ) && ( $an_bo eq '+' ) && ( $unan_bo eq '+' ) )
      ? '+'
      : ( ( $an_o eq '+' ) && ( $an_bo eq '+' ) && ( $unan_bo eq '-' ) ) ? '-'
      : ( ( $an_o eq '+' ) && ( $an_bo eq '-' ) && ( $unan_bo eq '+' ) ) ? '-'
      : ( ( $an_o eq '+' ) && ( $an_bo eq '-' ) && ( $unan_bo eq '-' ) ) ? '+'
      : ( ( $an_o eq '-' ) && ( $an_bo eq '+' ) && ( $unan_bo eq '+' ) ) ? '-'
      : ( ( $an_o eq '-' ) && ( $an_bo eq '+' ) && ( $unan_bo eq '-' ) ) ? '+'
      : ( ( $an_o eq '-' ) && ( $an_bo eq '-' ) && ( $unan_bo eq '+' ) ) ? '+'
      : ( ( $an_o eq '-' ) && ( $an_bo eq '-' ) && ( $unan_bo eq '-' ) ) ? '-'
      :                                                                    '+';
}

sub output_scf {
    my (
        $scf,            $gap_type, $chr,    $cM,
        $chr_loc_ref,    $length,   $orient, $comment,
        $scf_output_ref, $scf_bridges_ref
    ) = @_;
    my $separator =
      ( $gap_type eq 'recombination' )
      ? "\t--------------------------------------------------"
      : "";

    if ( !defined $scf_output_ref->{$scf} ) {
        if ( $gap_type ne 'ignore' ) {

            # Print gap
            print "$chr\t$chr_loc_ref->{pos}";
            $chr_loc_ref->{pos} += $gap_size{$gap_type} - 1;
            print
"\t$chr_loc_ref->{pos}\t$chr_loc_ref->{component}\tN\t$gap_size{$gap_type}\tfragment\tno\t# $gap_type$separator\n";
            $chr_loc_ref->{component}++;
            $chr_loc_ref->{pos}++;
        }

        # Print scaffold
        print "$chr\t$chr_loc_ref->{pos}";
        $chr_loc_ref->{pos} += $length - 1;
        print
"\t$chr_loc_ref->{pos}\t$chr_loc_ref->{component}\tD\t$scf\t1\t$length\t$orient";

        if ( !defined $comment ) {
            $comment = "# Added to $cM cM based on mate pairs";
        }
        print "\t$comment\n";

        $chr_loc_ref->{component}++;
        $chr_loc_ref->{pos}++;
        $scf_output_ref->{$scf}++;
    }

    return;
}
