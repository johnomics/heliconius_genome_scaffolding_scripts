#!/usr/bin/env perl

# convert_vcf_to_joinmap_heliconius_psti.pl

# This script is provided FOR INFORMATION ONLY and is not intended for use.
# It was written specifically for the Heliconius melpomene genome paper
# (Heliconius Genome Consortium, doi: 10.1038/nature11041) and has not
# been adapted for general use.

# Purpose: Generate markers suitable for linkage mapping based on
#          SNPs in a VCF file and list locations of markers on scaffolds
# Input : VCF file (v4.1) containing all individuals from an intercross
# Output: stdout: JoinMap-format markers for linkage mapping
#         stderr: scaffold positions and statistics for each marker

# Author: John Davey john.davey@ed.ac.uk
# Begun 07/07/11

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

my @offspring = (
    1,   2,   3,   4,   6,   7,   8,  12,  13,  19,  23,  24,
    25,  26,  28,  29,  30,  33,  35, 37,  39,  41,  57,  58,
    72,  73,  74,  75,  77,  78,  79, 100, 101, 103, 109, 110,
    111, 114, 115, 116, 117, 118, 119
);

my @females = (
    1,  4,  6,  7,  8,  19, 23, 24,  26,  30,
    33, 35, 57, 58, 72, 75, 77, 115, 116, 119
);
my @males = (
    2,  3,  12,  13,  25,  28,  29,  37,  39,  41,  73, 74,
    78, 79, 100, 101, 103, 109, 110, 111, 114, 117, 118
);

my $vcf_filename = "";
my $qthreshold   = 10;
my $mother_name  = "";
my $father_name  = "";
my $options_okay = GetOptions(
    'vcf=s'         => \$vcf_filename,
    'quality=i'     => \$qthreshold,
    'mother_name=s' => \$mother_name,
    'father_name=s' => \$father_name,
);
croak "No VCF file! Please specify -v $OS_ERROR\n" if ( $vcf_filename eq "" );
croak "No mother name! Please specify -m $OS_ERROR\n" if ( $mother_name eq "" );
croak "No father name! Please specify -m $OS_ERROR\n" if ( $father_name eq "" );

open my $vcf_file, '<', $vcf_filename
  or croak "Can't open $vcf_filename $OS_ERROR!\n";

my @sample_names;

my %marker_count;

my %sex;
map { $sex{$_} = 'f'; } @females;
map { $sex{$_} = 'm'; } @males;

my %joinmap;
my $bases_w_complete_calls_above_qthres = 0;
while ( my $vcf_line = <$vcf_file> ) {
    chomp $vcf_line;
    if ( $vcf_line =~ /^#CHROM/ ) {
        @sample_names = split /\t/, $vcf_line;
        for my $i ( 0 .. 8 ) { shift @sample_names; }
        map { s/PstI\.// } @sample_names;
    }
    next if ( $vcf_line =~ /^#/ );

    next if ( $vcf_line =~ /\.\/\./ );    # Discard bases with missing genotypes
    my %base;
    my @fields = split /\t/, $vcf_line;
    next if ( $fields[4] eq "." );        # If no alternate call
    my $skip_line = 0;
    for my $i ( 9 .. ( @sample_names + 8 ) ) {
        if ( $fields[$i] eq "./." ) { $skip_line++; last; }
        my @sample_fields = split /:/, $fields[$i];

        # Only works if GQ is the fourth field
        if ( $sample_fields[3] < $qthreshold ) { $skip_line++; last; }

        $base{ $sample_names[ $i - 9 ] }{gt} = $sample_fields[0];
        $base{ $sample_names[ $i - 9 ] }{dp} = $sample_fields[2];
        $base{ $sample_names[ $i - 9 ] }{gq} = $sample_fields[3];
    }

    next if ($skip_line);
    $bases_w_complete_calls_above_qthres++;
    my $f1pattern = "$base{$father_name}{gt} $base{$mother_name}{gt}";
    my $marker_type;
    my %f2_genotypes;
    if ( $f1pattern eq "0/0 0/1" ) {
        $marker_type         = "<lmxll>";
        $f2_genotypes{"0/0"} = "ll";
        $f2_genotypes{"0/1"} = "lm";
    }
    elsif ( $f1pattern eq "1/1 0/1" ) {
        $marker_type = "<lmxll>";
        $f2_genotypes{"0/1"} = "lm";
        $f2_genotypes{"1/1"} = "ll";
    }
    elsif ( $f1pattern eq "0/1 0/0" ) {
        $marker_type         = "<nnxnp>";
        $f2_genotypes{"0/0"} = "nn";
        $f2_genotypes{"0/1"} = "np";
        $f2_genotypes{"1/1"} = "p-";        # Sex-linked only
    }
    elsif ( $f1pattern eq "0/1 1/1" ) {
        $marker_type = "<nnxnp>";

        $f2_genotypes{"0/0"} = "p-";        # Sex-linked only
        $f2_genotypes{"0/1"} = "np";
        $f2_genotypes{"1/1"} = "nn";
    }
    elsif ( $f1pattern eq "0/1 0/1" ) {
        $marker_type         = "<hkxhk>";
        $f2_genotypes{"0/0"} = "hh";
        $f2_genotypes{"0/1"} = "hk";
        $f2_genotypes{"1/1"} = "kk";
    }
    elsif ( $f1pattern eq "0/0 1/1" ) {
        $marker_type         = "sex_aa:b-";
        $f2_genotypes{"0/1"} = "lm";
        $f2_genotypes{"0/0"} = "ll";
    }
    elsif ( $f1pattern eq "1/1 0/0" ) {
        $marker_type         = "sex_aa:b-";
        $f2_genotypes{"0/1"} = "lm";
        $f2_genotypes{"1/1"} = "ll";
    }
    else {
        next;
    }
    my $scaffold                 = $fields[0];
    my $joinmap_marker           = "";
    my $sex_linked_nnxnp_males   = 0;
    my $sex_linked_nnxnp_females = 0;
    foreach my $offspring (@offspring) {
        my $gt = $f2_genotypes{ $base{$offspring}{gt} };
        if (   ( defined $base{$offspring}{gt} )
            && ( defined $gt ) )
        {
            if ( $marker_type eq "<nnxnp>" ) {
                if ( ( $sex{$offspring} eq 'm' ) && ( $gt eq 'p-' ) ) {
                    $skip_line++;
                    last;
                }
                if (   ( $sex{$offspring} eq 'm' )
                    && ( ( $gt eq 'nn' ) || ( $gt eq 'np' ) ) )
                {
                    $sex_linked_nnxnp_males++;
                }
                if (( $sex{$offspring} eq 'f' )
                  && ( ( $gt eq 'nn' ) || ( $gt eq 'p-' ) ) ){
                      $sex_linked_nnxnp_females++;
                  }

                  $joinmap_marker .= "$gt ";
            }
            else {
                  $joinmap_marker .= "$gt ";
            }
        }
        else {
              $skip_line++;
              last;
        }
    }
    next if ($skip_line);

    if (     ( $sex_linked_nnxnp_males == @males )
          && ( $sex_linked_nnxnp_females == @females ) )
      {
          $marker_type = "sex_ab:a-";
          my @gts = split / /, $joinmap_marker;
          $joinmap_marker = "";
          for my $i ( 0 .. $#gts ) {
              if ( ( $sex{ $offspring[$i] } eq 'f' ) && ( $gts[$i] eq 'nn' ) ) {
                  $gts[$i] = 'n-';
              }
              $joinmap_marker .= "$gts[$i] ";
          }
    }
    else {
        next if ($joinmap_marker =~ /p\-/);
    }

    my $male_depth = 0;
    my $male_qual  = 0;

    my $average_depth = 0;
    my $average_qual  = 0;

    foreach my $male (@males) {
          $average_depth += $base{$male}{dp};
          $average_qual  += $base{$male}{gq};
          $male_depth    += $base{$male}{dp};
          $male_qual     += $base{$male}{gq};
    }

    my $female_depth = 0;
    my $female_qual  = 0;
    foreach my $female (@females) {
          $average_depth += $base{$female}{dp};
          $average_qual  += $base{$female}{gq};
          $female_depth  += $base{$female}{dp};
          $female_qual   += $base{$female}{gq};
    }

    my $pos_dp = 0;
    my $pos_mq = 0;
    if ( $fields[7] =~ /DP=(\d+?);/ ) {
          $pos_dp = $1;
    }
    if ( $fields[7] =~ /MQ=(.+?);/ ) {
          $pos_mq = $1;
    }
    chop $joinmap_marker;
    $average_depth /= @offspring;
    $average_qual  /= @offspring;
    $female_depth  /= @females;
    $female_qual   /= @females;
    $male_depth    /= @males;
    $male_qual     /= @males;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{dp} =
      int $average_depth;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{gq} =
      int $average_qual;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{pos_dp} =
      $pos_dp;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{pos_mq} =
      $pos_mq;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{male_dp} =
      int $male_depth;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{male_gq} =
      int $male_qual;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{female_dp} =
      int $female_depth;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{female_gq} =
      int $female_qual;
}

close $vcf_file;

my $marker_num = 0;
foreach my $type ( sort keys %joinmap ) {
      my $scf_marker = 0;
      foreach my $marker ( keys %{ $joinmap{$type} } ) {

          # Output markers appearing on more than one scaffold
          # or appearing more than once on one scaffold
          my $check_marker_num = 0;
          my $num_bases        = scalar keys %{ $joinmap{$type}{$marker} };
          next if ( $num_bases <= 1 );
          $marker_num++;
          my $sex_indicator = "";
          my $output_type   = $type;
          if ( $type eq "sex_aa:b-" ) {
              $output_type   = "<lmxll>";
              $sex_indicator = "_sex_aa:b-";
          }
          elsif ( $type eq "sex_ab:a-" ) {
              $output_type   = "<nnxnp>";
              $sex_indicator = "_sex_ab:a-";
          }
          print "$marker_num$sex_indicator:$num_bases\t$output_type\t$marker\n";
          foreach my $base ( sort keys %{ $joinmap{$type}{$marker} } ) {
              print STDERR
"$marker_num$sex_indicator\t$base\t$output_type\t$marker\t$joinmap{$type}{$marker}{$base}{dp}\t$joinmap{$type}{$marker}{$base}{gq}\t$joinmap{$type}{$marker}{$base}{pos_dp}\t$joinmap{$type}{$marker}{$base}{pos_mq}\t$joinmap{$type}{$marker}{$base}{male_dp}\t$joinmap{$type}{$marker}{$base}{male_gq}\t$joinmap{$type}{$marker}{$base}{female_dp}\t$joinmap{$type}{$marker}{$base}{female_gq}\n";
          }
      }
}
