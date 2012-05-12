#!/usr/bin/env perl

# validate_markers_for_mapping.pl

# This script is provided FOR INFORMATION ONLY and is not intended for use.
# It was written specifically for the Heliconius melpomene genome paper
# (Heliconius Genome Consortium, doi: 10.1038/nature11041) and has not
# been adapted for general use.

# Purpose: Load marker positions and validate markers:
#          - error-correct female-specific chromosome prints
#          - convert hkxhk markers to nnxnp markers wherever possible
#          - output remaining nnxnp markers
#          - output remaining hkxhk markers

# Input : marker position list from convert_vcf_to_joinmap_heliconius.pl
# Output: valid lmxll markers
#         valid nnxnp markers (including converted hkxhk markers)
#         remaining hkxhk markers
#         remaining nnxnp markers

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

# Autoflush output so reporting on progress works
$| = 1;

my $marker_filename = "";

my $options_okay = GetOptions(
    'markers=s' => \$marker_filename,
);
croak "No markers file! Please specify -m $OS_ERROR\n"
  if ( $marker_filename eq "" );


# Load markers
open my $marker_file, '<', $marker_filename
  or croak "Can't open $marker_filename $OS_ERROR!\n";

my %gt;
my $nind = 0;
while ( my $marker_line = <$marker_file> ) {
    chomp $marker_line;
    my ( $marker, $scfpos, $type, $genotypes ) = split /\t/, $marker_line;
    $gt{$type}{$marker}{gt} = $genotypes;
    if ($nind == 0) {
        my @gts = split / /, $genotypes;
        $nind = @gts;
    }
    $gt{$type}{$marker}{count}++;
}
close $marker_file;

# Collapse mirror and error female-specific patterns
my %ll_lookup;
my %ll_merged;

my @sorted_markers =
  reverse sort { $gt{"<lmxll>"}{$a}{count} <=> $gt{"<lmxll>"}{$b}{count} }
  keys %{ $gt{"<lmxll>"} };

foreach my $marker (@sorted_markers) {
    my @gts     = split / /, $gt{"<lmxll>"}{$marker}{gt};
    my $pattern = "*$marker\_$gt{\"<lmxll>\"}{$marker}{count}\_<lmxll>\t";
    my $mirror  = "*$marker\_$gt{\"<lmxll>\"}{$marker}{count}\_<lmxll>_M\t";
    foreach my $gt (@gts) {
        $pattern .= ( $gt eq "nn" ) ? "A " : "H ";
        $mirror  .= ( $gt eq "nn" ) ? "H " : "A ";
    }
    chop $pattern;
    chop $mirror;

}
my $unfiltered_lmxll = scalar keys %{ $gt{"<lmxll>"} };

foreach my $marker (@sorted_markers) {
    if ( defined $gt{"<lmxll>"}{$marker} ) {
        $ll_lookup{$marker} = $marker;
        $gt{"<lmxll>"}{$marker}{lookup}{$marker}++;
        foreach my $marker2 (@sorted_markers) {

            next if ( !defined $gt{"<lmxll>"}{$marker2} );

            # Ignore marker2 if already processed
            next
              if ( $gt{"<lmxll>"}{$marker2}{count} >
                $gt{"<lmxll>"}{$marker}{count} );

            next if ( $marker2 eq $marker );
            if (
                is_match(
                    $gt{"<lmxll>"}{$marker}{gt},
                    $gt{"<lmxll>"}{$marker2}{gt}
                )
              )
            {
                print "$marker matches $marker2\n";
                $ll_lookup{$marker2} = $marker;
                $ll_merged{$marker2}{count} = $gt{"<lmxll>"}{$marker2}{count};
                $ll_merged{$marker2}{gt} = $gt{"<lmxll>"}{$marker2}{gt};
                $gt{"<lmxll>"}{$marker}{lookup}{$marker2}++;
                $gt{"<lmxll>"}{$marker}{count} +=
                  $gt{"<lmxll>"}{$marker2}{count};
                delete $gt{"<lmxll>"}{$marker2};
            }
        }
    }
}
print STDERR "$unfiltered_lmxll lmxll markers collapsed to ";
print STDERR scalar keys %{ $gt{"<lmxll>"} };
print STDERR "\n";

sub is_match {
    my ( $marker1, $marker2 ) = @_;
    my @gt1 = split /\s/, $marker1;
    my @gt2 = split /\s/, $marker2;
    return 0 if ( @gt1 != @gt2 );

    my $mirror_mismatch = 0;
    my $same_mismatch   = 0;
    for my $i ( 0 .. $#gt1 ) {
        $same_mismatch++ if ( $gt1[$i] ne $gt2[$i] );
        my $gtm = $gt2[$i] eq "lm" ? "ll" : "lm";
        $mirror_mismatch++ if ( $gt1[$i] ne $gtm );
        return 0 if ( ( $mirror_mismatch > 3 ) && ( $same_mismatch > 3 ) );
    }
    return 1;
}

my $unfiltered_hkxhk = scalar keys %{ $gt{"<hkxhk>"} };
foreach my $hh_marker ( keys %{ $gt{"<hkxhk>"} } ) {
    my @gts = split / /, $gt{"<hkxhk>"}{$hh_marker}{gt};
    my %gts;
    $gts{hh} = 0;
    $gts{hk} = 0;
    $gts{kk} = 0;
    map { $gts{$_}++; } @gts;
    if (
        ( chi_sq( $gts{hh}, $gts{hk}, $gts{kk} ) > 3.84 )
    or ( ( $gts{hh} == 0 ) or ( $gts{kk} == 0 ) )
      )
    {
        delete $gt{"<hkxhk>"}{$hh_marker};
        next;
    }
}

sub chi_sq {
    my ( $aa, $ab, $bb ) = @_;
    my $p      = ( 2 * $aa + $ab ) / ( 2 * ( $aa + $ab + $bb ) );
    my $q      = 1 - $p;
    my $n      = $aa + $ab + $bb;
    my $exp_aa = $p**2 * $n;
    my $exp_ab = 2 * $p * $q * $n;
    my $exp_bb = $q**2 * $n;

    my $chi_sq =
      ( $aa - $exp_aa )**2 / $exp_aa +
      ( $ab - $exp_ab )**2 / $exp_ab +
      ( $bb - $exp_bb )**2 / $exp_bb;
    return $chi_sq;
}

print STDERR "$unfiltered_hkxhk hkxhk markers collapsed to ";
print STDERR scalar keys %{ $gt{"<hkxhk>"} };
print STDERR "\n";

my %ll_links;
my @ll_out_lines;
my %ll_chrom_marker;
foreach my $ll_marker ( sort { (split /_/,$a)[0] <=> (split /_/,$b)[0] } keys %{ $gt{"<lmxll>"} } ) {
    my $ll_hhkk = $gt{"<lmxll>"}{$ll_marker}{gt};
    my $ll_kkhh = $gt{"<lmxll>"}{$ll_marker}{gt};

    $ll_hhkk =~ s/ll/hh/g;
    $ll_hhkk =~ s/lm/kk/g;
    $ll_kkhh =~ s/ll/kk/g;
    $ll_kkhh =~ s/lm/hh/g;

    my $hh_found = 0;
    foreach my $hh_marker ( keys %{ $gt{"<hkxhk>"} } ) {

        if (   ( hk_match( $ll_hhkk, $gt{"<hkxhk>"}{$hh_marker}{gt} ) )
            or ( hk_match( $ll_kkhh, $gt{"<hkxhk>"}{$hh_marker}{gt} ) ) )
        {
            $hh_found++;
            $ll_links{$hh_marker}{$ll_marker}++;
        }
    }

    # If no matching hkxhk markers for this lmxll marker, delete it
    if ( $hh_found == 0 ) {
        foreach
          my $ll_lookup_marker ( keys %{ $gt{"<lmxll>"}{$ll_marker}{lookup} } )
        {
            delete $ll_lookup{$ll_lookup_marker};
        }
        delete $gt{"<lmxll>"}{$ll_marker};
    }
    else {
        push @ll_out_lines, "$ll_marker\_$gt{\"<lmxll>\"}{$ll_marker}{count}\_$hh_found\_$ll_lookup{$ll_marker}\t<lmxll>\t$gt{\"<lmxll>\"}{$ll_marker}{gt}\n";
        $ll_chrom_marker{$ll_marker}++;
    }
}

print "name=$marker_filename\_lmxll\n";
print "popt=CP\n";
print "nloc=" . (scalar @ll_out_lines) . "\n";
print "nind=$nind\n";
print @ll_out_lines;

foreach my $ll_marker (keys %ll_lookup) {
    if (!defined ($ll_chrom_marker{$ll_marker})) {
       print "$ll_marker\_$ll_merged{$ll_marker}{count}\_0\_$ll_lookup{$ll_marker}\t<lmxll>\t$ll_merged{$ll_marker}{gt}\n";
    }
}

sub hk_match {
    my ( $full, $partial ) = @_;
    my @full    = split / /, $full;
    my @partial = split / /, $partial;

    return 0 if ( @full != @partial );

    foreach my $i ( 0 .. $#full ) {
        next if ( $partial[$i] eq "hk" );
        return 0 if ( $full[$i] ne $partial[$i] );
    }
    return 1;
}

print STDERR scalar keys %{ $gt{"<lmxll>"} };
print STDERR " lmxll markers with links to ";
print STDERR scalar keys %ll_links;
print STDERR " hkxhk markers\n";


my %hh_nn_markers;
my @hh_no_ll_markers;
my @hh_nn_markers;
foreach my $hh_marker ( sort { $a <=> $b } keys %{ $gt{"<hkxhk>"} } ) {

    if ( !defined( $ll_links{$hh_marker} ) ) {
        push @hh_no_ll_markers,
"$hh_marker\_$gt{\"<hkxhk>\"}{$hh_marker}{count}\t<hkxhk>\t$gt{\"<hkxhk>\"}{$hh_marker}{gt}\n";
        next;
    }
    foreach my $ll_marker ( keys %{ $ll_links{$hh_marker} } ) {
        if ( !( defined( $gt{"<lmxll>"}{$ll_marker} ) ) ) {
            delete $ll_links{$hh_marker}{$ll_marker};
        }
        my $nn_marker = "";
        my @ll_gts    = split / /, $gt{"<lmxll>"}{$ll_marker}{gt};
        my @hh_gts    = split / /, $gt{"<hkxhk>"}{$hh_marker}{gt};
        next if ( @ll_gts != @hh_gts );
        my %hh_ll_pair_check;
        foreach my $i ( 0 .. $#ll_gts ) {
            $hh_ll_pair_check{ $hh_gts[$i] }{ $ll_gts[$i] }++;
        }
        my $hh_l_match = (
            reverse sort {
                $hh_ll_pair_check{"hh"}{$a} <=> $hh_ll_pair_check{"hh"}{$b}
              } keys %{ $hh_ll_pair_check{"hh"} }
        )[0];
        my $kk_l_match = (
            reverse sort {
                $hh_ll_pair_check{"kk"}{$a} <=> $hh_ll_pair_check{"kk"}{$b}
              } keys %{ $hh_ll_pair_check{"kk"} }
        )[0];
        next if ( $hh_l_match eq $kk_l_match );

        foreach my $i ( 0 .. $#ll_gts ) {
            $nn_marker .=
              ( ( $hh_gts[$i] eq "hh" ) && ( $ll_gts[$i] eq $hh_l_match ) )
              ? "nn "
              : ( ( $hh_gts[$i] eq "kk" ) && ( $ll_gts[$i] eq $kk_l_match ) )
              ? "np "
              : ( ( $hh_gts[$i] eq "hk" ) && ( $ll_gts[$i] eq $hh_l_match ) )
              ? "np "
              : ( ( $hh_gts[$i] eq "hk" ) && ( $ll_gts[$i] eq $kk_l_match ) )
              ? "nn "
              : "-- ";
        }
        chop $nn_marker;    #Remove trailing space
        $hh_nn_markers{$nn_marker}{$ll_marker}{$hh_marker}++;
        push @hh_nn_markers, 
"$hh_marker\_$gt{\"<hkxhk>\"}{$hh_marker}{count}\_$ll_marker\t<nnxnp>\t$nn_marker\n";
    }
}
print "name=$marker_filename\_hkxhk_converted_to_nnxnp\n";
print "popt=CP\n";
print "nloc=" . (scalar @hh_nn_markers) . "\n";
print "nind=$nind\n";
print @hh_nn_markers;

print "name=$marker_filename\_hkxhk_unlinked_to_lmxll\n";
print "popt=CP\n";
print "nloc=" . (scalar @hh_no_ll_markers) . "\n";
print "nind=$nind\n";
print @hh_no_ll_markers;


my @nn_no_hk_markers;
foreach my $nn_marker (
    sort { $gt{"<nnxnp>"}{$a} <=> $gt{"<nnxnp>"}{$b} }
    keys %{ $gt{"<nnxnp>"} }
  )
{
    my $nn_unmatched = 1;
    foreach my $hh_nn_marker ( keys %hh_nn_markers ) {
        next
          if ( keys %{ $hh_nn_markers{$hh_nn_marker} } > 1 )
          ;    # Ignore if nn marker matches more than one ll marker
        if ( $gt{"<nnxnp>"}{$nn_marker}{gt} eq $hh_nn_marker ) {
            $ll_links{$nn_marker}
              { ( keys %{ $hh_nn_markers{$hh_nn_marker} } )[0] }++;
            $nn_unmatched = 0;
        }
    }
    if ($nn_unmatched) {
        push @nn_no_hk_markers,
"$nn_marker\_$gt{\"<nnxnp>\"}{$nn_marker}{count}\t<nnxnp>\t$gt{\"<nnxnp>\"}{$nn_marker}{gt}\n";
    }
}

print "name=$marker_filename\_nnxnp_unlinked_to_hkxhk\n";
print "popt=CP\n";
print "nloc=" . (scalar @nn_no_hk_markers) . "\n";
print "nind=$nind\n";
print @nn_no_hk_markers;