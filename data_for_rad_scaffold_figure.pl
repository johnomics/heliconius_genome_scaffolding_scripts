#!/usr/bin/env perl

# data_for_rad_scaffold_figure.pl

# This script is provided FOR INFORMATION ONLY and is not intended for use.
# It was written specifically for the Heliconius melpomene genome paper
# (Heliconius Genome Consortium, doi: 10.1038/nature11041) and has not
# been adapted for general use.

# Purpose: Generate input data for draw_rad_scaffold_figure.R,
#          used to produce Supplementary Figure S4.6.1

# Input  : chromosome to output
#          chromosome AGP file
#          scaffold AGP file
#          full scaffolding information from scaffold_heliconius_genome.pl
# Output : scaffold start and end positions
#          marker positions
#          cM positions

# Author: John Davey john.davey@ed.ac.uk
# Begun 23/10/2011

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

my $chromosome_filename  = "";
my $scaffold_filename    = "";
my $full_output_filename = "";
my $chromosome           = 18;

my $options_okay = GetOptions(
    'agp_chrom=s'   => \$chromosome_filename,
    'scaffolds=s'   => \$scaffold_filename,
    'full_output=s' => \$full_output_filename,
    'chromosome=s'  => \$chromosome,
);

croak
"\nUsage: perl make_hgp_figure2_input.pl -a chromosome_agp_file -s scaffold_agp_file -f full_output_file -c chromosome_num\n"
  if !$options_okay;

croak "No chromosome filename! $OS_ERROR\n"  if ( $chromosome_filename  eq "" );
croak "No scaffold filename! $OS_ERROR\n"    if ( $scaffold_filename    eq "" );
croak "No full output filename! $OS_ERROR\n" if ( $full_output_filename eq "" );

my %chr_scfs;
my $chr_length = 0;
open my $chromosome_file, '<', $chromosome_filename
  or croak "Can't open $chromosome_filename: $OS_ERROR!\n";
while ( my $chromosome_line = <$chromosome_file> ) {
    if ( $chromosome_line =~ /^chr$chromosome\t(.+)scf/ ) {
        my @f = split /\t/, $chromosome_line;
        my $scf = $f[5];
        $chr_scfs{$scf}{start} = $f[1];
        $chr_scfs{$scf}{end}   = $f[2];
        $chr_scfs{$scf}{dir}   = $f[8];
        $chr_length            = $chr_scfs{$scf}{end};
    }
}
close $chromosome_file;

open my $scf_out_file, '>', "chr$chromosome\.scf.start.end.txt"
  or croak "Can't open chr$chromosome\.scf.start.end.txt: $OS_ERROR!\n";
print $scf_out_file "Scaffold\tStart\tEnd\tDir\n";
foreach my $scf ( sort { $chr_scfs{$a}{start} <=> $chr_scfs{$b}{start} }
    keys %chr_scfs )
{
    print $scf_out_file
"$scf\t$chr_scfs{$scf}{start}\t$chr_scfs{$scf}{end}\t$chr_scfs{$scf}{dir}\n";
}
close $scf_out_file;

open my $scaffold_file, '<', $scaffold_filename
  or croak "Can't open $scaffold_filename: $OS_ERROR!\n";
my %scf_ctgs;
while ( my $scaffold_line = <$scaffold_file> ) {
    chomp $scaffold_line;
    my @f = split /\t/, $scaffold_line;
    my $scf = $f[0];
    if ( defined $chr_scfs{$scf} ) {
        next if ( $f[4] eq 'N' );
        my $ctg = $f[3];
        $scf_ctgs{$scf}{$ctg}{start} = $f[1];
        $scf_ctgs{$scf}{$ctg}{end}   = $f[2];
        if ( $chr_scfs{$scf}{dir} eq '+' ) {
            $scf_ctgs{$scf}{$ctg}{chrstart} =
              $chr_scfs{$scf}{start} + $f[1] - 1;
            $scf_ctgs{$scf}{$ctg}{chrend} = $chr_scfs{$scf}{start} + $f[2] - 1;
        }
        else {
            $scf_ctgs{$scf}{$ctg}{chrstart} = $chr_scfs{$scf}{end} - $f[1] + 1;
            $scf_ctgs{$scf}{$ctg}{chrend}   = $chr_scfs{$scf}{end} - $f[2] + 1;
        }
    }
}

my %marker;

open my $full_output_file, '<', $full_output_filename
  or croak "Can't open $full_output_filename: $OS_ERROR!\n";
while ( my $full_output_line = <$full_output_file> ) {
    chomp $full_output_line;
    next if ( $full_output_line =~ /HMEL/ );
    my @f = split /\t/, $full_output_line;
    if ( ( defined $f[0] ) && ( defined $chr_scfs{ $f[0] } ) ) {
        my $pos    = $f[2];
        my $scf    = $f[0];
        my $chrpos = 0;
        if ( $chr_scfs{$scf}{dir} eq '+' ) {
            $chrpos = $chr_scfs{$scf}{start} + $pos - 1;
        }
        else {
            $chrpos = $chr_scfs{$scf}{end} - $pos + 1;
        }
        $marker{$chrpos}{scf} = $f[0];
        $marker{$chrpos}{ctg} = $f[1];
        $marker{$chrpos}{cm}  = $f[4];
        $marker{$chrpos}{pos} = $pos;
    }
}

open my $marker_out_file, '>', "chr$chromosome.scf.marker.pos.cm.txt"
  or croak "Can't open chr$chromosome.scf.marker.pos.cm.txt: $OS_ERROR\n";
print $marker_out_file "Scaffold\tContig\tChrom Position\tScaff Position\tcM\n";
foreach my $chrpos ( sort { $a <=> $b } keys %marker ) {
        print $marker_out_file
          "$marker{$chrpos}{scf}\t$marker{$chrpos}{ctg}\t$chrpos\t$marker{$chrpos}{pos}\t";
        print $marker_out_file $marker{$chrpos}{cm} eq "-"
          ? "-1"
          : $marker{$chrpos}{cm};
        print $marker_out_file "\n";
}
close $marker_out_file;

my $last_cm  = "0.000";
my $last_pos = 0;
my $last_scf = "";
my $last_ctg = "";
my %cmpos;
$cmpos{$last_cm}{start} = 1;

foreach my $pos ( sort { $a <=> $b } keys %marker ) {

    if ( $marker{$pos}{cm} ne '-' ) {
        if ( $marker{$pos}{cm} ne $last_cm ) {
            my $scf        = $marker{$pos}{scf};
            my $ctg        = $marker{$pos}{ctg};
            my $breakpoint = 0;
            if ( ( $scf eq $last_scf ) and ( $ctg eq $last_ctg ) ) {
                $breakpoint = int( $last_pos + ( $pos - $last_pos ) / 2 );
                $cmpos{$last_cm}{end} = $breakpoint;
                $cmpos{ $marker{$pos}{cm} }{start} = $breakpoint + 1;
            }
            else {
                $cmpos{$last_cm}{end} =
                    $chr_scfs{$last_scf}{dir} eq '+'
                  ? $scf_ctgs{$last_scf}{$last_ctg}{chrend}
                  : $scf_ctgs{$last_scf}{$last_ctg}{chrstart};
                $cmpos{ $marker{$pos}{cm} }{start} =
                    $chr_scfs{$scf}{dir} eq '+'
                  ? $scf_ctgs{$scf}{$ctg}{chrstart}
                  : $scf_ctgs{$scf}{$ctg}{chrend};
            }
        }
        $last_cm  = $marker{$pos}{cm};
        $last_scf = $marker{$pos}{scf};
        $last_ctg = $marker{$pos}{ctg};

        $last_pos = $pos;
    }
}
$cmpos{$last_cm}{end} = $chr_length;

open my $cmpos_out_file, '>', "chr$chromosome\.cm.pos.txt"
  or croak "Can't open chr$chromosome\.cm.pos.txt: $OS_ERROR!\n";
print $cmpos_out_file "cM\tStart\tEnd\n";
foreach my $cm ( sort { $a <=> $b } keys %cmpos ) {
    print $cmpos_out_file "$cm\t$cmpos{$cm}{start}\t$cmpos{$cm}{end}\n";
}
close $cmpos_out_file;
