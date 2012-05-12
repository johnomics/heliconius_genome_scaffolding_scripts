#!/usr/bin/env perl

# generate_bambus2_agps.pl
#
# This script is provided FOR INFORMATION ONLY and is not intended for use.
# It was written specifically for the Heliconius melpomene genome paper
# (Heliconius Genome Consortium, doi: 10.1038/nature11041) and has not
# been adapted for general use.

# Purpose: Run Bambus2 for each chromosome to identify local matepair bridges

# Input  : chromosome AGP from scaffold_heliconius_genome.pl
#          revised scaffold AGP from scaffold_heliconius_genome.pl
#          tab-delimited list of matepair alignments, eg
#              MatePairName\tScaffoldName\tStartPos\tEndPos\tDir
#            where MatePairName ends with 'a' or 'b'
#            and Dir is 'f' or 'r' for forward / reverse
#          tab-delimited list of genome scaffold lengths, eg
#            scf00001\t635048
#          genome contigs in FASTA format

# Output : folder of Bambus input files and Bambus-generated AGP files
#          info files recording usage of matepair reads for scaffolding

# Requires Bambus2 (written for AMOS 3.0.1)

# Author: John Davey john.davey@ed.ac.uk

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
use Graph;

# Autoflush output so reporting on progress works
$| = 1;

my $chr_agp_in = "";
my $scf_agp_in = "";
my $matepair_filename = "";
my $lengths_filename  = "";
my $contigs_filename  = "";

my $dir_name = "bambus_input";
my $debug;
my $verbose;


my $options_okay = GetOptions(
    'chr_in=s'    => \$chr_agp_in,
    'scf_in=s'    => \$scf_agp_in,
    'matepairs=s' => \$matepair_filename,
    'lengths=s'   => \$lengths_filename,
    'contigs=s'   => \$contigs_filename,
    'dir=s'       => \$dir_name,
    'debug'       => \$debug,
    'verbose'     => \$verbose,
);

croak "No input chromosome AGP file! Please specify --chr_in $OS_ERROR\n"
  if ( $chr_agp_in eq "" );
croak "No input scaffold AGP file! Please specify --scf_in $OS_ERROR\n"
  if ( $scf_agp_in eq "" );
croak "No matepairs file! Please specify --matepairs $OS_ERROR\n"
  if ( $matepair_filename eq "" );
croak "No lengths file! Please specify --lengths $OS_ERROR\n"
  if ( $lengths_filename eq "" );
croak "No contigs file! Please specify --contigs $OS_ERROR\n"
  if ( $contigs_filename eq "" );
croak "No output directory! Please specify --dir $OS_ERROR\n"
  if ( $dir_name eq "" );

print STDERR "Loading contigs...\n";
my %contigs;
open my $contigs_file, '<', $contigs_filename
  or croak "Can't open $contigs_filename $OS_ERROR!\n";
while ( my $contigs_line = <$contigs_file> ) {
    chomp $contigs_line;
    if ( $contigs_line =~ /^>(.+)/ ) {
        my $ctg_name = $1;
        my $ctg_seq  = <$contigs_file>;
        chomp $ctg_seq;
        $contigs{$ctg_name} = $ctg_seq;
    }
}

print STDERR "Loading scaffold AGP...\n";

# Define mate pair libraries
my %library;

$library{'3kb'}{min} = 700;
$library{'3kb'}{max} = 4342;
$library{'8kb'}{min} = 1665;
$library{'8kb'}{max} = 8331;

my %scf_contigs;
my %contigs_scf;
my %scf_not_on_map;
open my $scf_agp, '<', $scf_agp_in
  or croak "Can't open $scf_agp_in $OS_ERROR!\n";
while ( my $scf_agp_line = <$scf_agp> ) {
    chomp $scf_agp_line;
    my @agp_fields = split /\t/, $scf_agp_line;
    my $scf = $agp_fields[0];
    if ( $agp_fields[4] eq 'N' ) {
        my $ctg_id = $agp_fields[3] - 1;
        next if ( $ctg_id < 1 );
        if ( !defined $scf_contigs{$scf} ) {
            print STDERR "No scf for gap $scf $ctg_id!\n";
        }
        if ( !defined $scf_contigs{$scf}{$ctg_id} ) {
            print STDERR "No ctg id for gap $scf $ctg_id!\n";
        }
        $scf_contigs{$scf}{$ctg_id}{gap} = $agp_fields[5];
        next;
    }
    my $ctg_id   = $agp_fields[3];
    my $ctg_name = $agp_fields[5];
    $scf_contigs{$scf}{$ctg_id}{name}   = $ctg_name;
    $scf_contigs{$scf}{$ctg_id}{orient} = $agp_fields[8];
    $contigs_scf{$ctg_name}{scf}        = $scf;
    $contigs_scf{$ctg_name}{scf_start}  = $agp_fields[1];
    $contigs_scf{$ctg_name}{scf_end}    = $agp_fields[2];
    $contigs_scf{$ctg_name}{scf_orient} = $agp_fields[8];
    $contigs_scf{$ctg_name}{ctg_start}  = $agp_fields[6];
    $contigs_scf{$ctg_name}{ctg_end}    = $agp_fields[7];
    $scf_not_on_map{$scf}++;
}
close $scf_agp;

print STDERR "Loading scaffold lengths...\n";
my %scf_lengths;
open my $lengths_file, '<', $lengths_filename
  or croak "Can't open $lengths_filename $OS_ERROR!\n";
while ( my $lengths_line = <$lengths_file> ) {
    chomp $lengths_line;
    my ( $scf, $length ) = split /\t/, $lengths_line;
    $scf_lengths{$scf} = $length;
}

print STDERR "Loading matepairs...\n";
my %mates;
my %ctg_bridges;

open my $matepairs, "<", $matepair_filename
  or croak "Can't open $matepair_filename $OS_ERROR\n";
my $mate_count;
while ( my $matepair_line = <$matepairs> ) {
    chomp $matepair_line;
    my ( $matepair, $contig, $startpos, $endpos, $dir ) = split /\t/,
      $matepair_line;
    my $ctg_name = "ctg" . $contig;
    next if ( !defined( $contigs_scf{$ctg_name}{scf} ) );
    $mate_count++;
    if ( $mate_count % 10000 == 0 )  { print STDERR '.'; }
    if ( $mate_count % 100000 == 0 ) { print STDERR "$mate_count\n"; }

    my $readname = substr $matepair, 0, -1;
    my $end = substr $matepair, -1;
    my $library = $readname =~ /^F/ ? "3kb" : "8kb";
    $mates{$readname}{library} = $library;
    $mates{$readname}{ends}{$end}{ctg} = $ctg_name;

    if ( $contigs_scf{$ctg_name}{scf_orient} eq '+' ) {
        $mates{$readname}{ends}{$end}{scf_orient} =
          ( $dir eq 'f' ) ? '' : ( $dir eq 'r' ) ? 'RC' : '';
        $mates{$readname}{ends}{$end}{scf_start} =
          $contigs_scf{$ctg_name}{scf_start} + $startpos - 1;
        $mates{$readname}{ends}{$end}{scf_end} =
          $contigs_scf{$ctg_name}{scf_start} + $endpos - 1;
    }

    if ( $contigs_scf{$ctg_name}{scf_orient} eq '-' ) {
        $mates{$readname}{ends}{$end}{scf_orient} =
          ( $dir eq 'f' ) ? 'RC' : ( $dir eq 'r' ) ? '' : 'RC';
        $mates{$readname}{ends}{$end}{scf_start} =
          $contigs_scf{$ctg_name}{scf_end} - $endpos + 1;
        $mates{$readname}{ends}{$end}{scf_end} =
          $contigs_scf{$ctg_name}{scf_end} - $startpos + 1;
    }
    $ctg_bridges{$ctg_name}{$readname} = $end;
}
close $matepairs;

print STDERR "\nFound " . ( scalar keys %mates ) . " mate pairs\n";

print STDERR "Loading chromosome AGP...\n";
open my $chr_agp, '<', $chr_agp_in
  or croak "Can't open $chr_agp_in $OS_ERROR!\n";

my %chr_scfs;
my %chr_cM_bounds;
while ( my $chr_agp_line = <$chr_agp> ) {
    chomp $chr_agp_line;
    my @agp_fields = split /\t/, $chr_agp_line;
    next if ( $agp_fields[4] eq 'N' );
    my $chr           = $agp_fields[0];
    my $scf           = $agp_fields[5];
    my $scf_start_pos = $agp_fields[1];
    my $scf_end_pos   = $agp_fields[2];
    $scf_lengths{$scf} = $agp_fields[7];
    delete $scf_not_on_map{$scf};

    my @cMs;
    if ( $chr =~ /unmapped/ ) {

        $chr = ( split /_/, $chr )[0];    # Strip _unmapped
        push @cMs, "-1";
    }
    else {
        my @loci = split /\|/, $agp_fields[9];
        shift @loci;                      # Remove '# Contigs :'
        foreach my $locus (@loci) {
            if ( $locus =~ /\ (.+)\ cM:\ (.+)\ / ) {
                push @cMs, $1;
            }
        }
    }
    foreach my $cM (@cMs) {
        $chr_scfs{$chr}{$scf}{$cM}++;
        $chr_cM_bounds{$chr}{$cM}{$scf_start_pos} = $scf;
        $chr_cM_bounds{$chr}{$cM}{$scf_end_pos}   = $scf;
    }
}

close $chr_agp;

print STDERR "Processing chromosomes...\n";
my %all_reads;
my %same_scf_reads;
my %spanning_reads;
my %missing_end_reads;
my %off_chrom_reads;
my %off_map_reads;
my %distant_reads;
my %read_seen;

mkdir $dir_name;

foreach my $chr ( sort keys %chr_scfs ) {

    print STDERR "$chr\tGenerate input";

    open my $linkage_xml, '>', "$dir_name/$chr.linkage.xml"
      or croak "Can't open linkage evidence XML file for $chr $OS_ERROR\n";
    print $linkage_xml "<EVIDENCE>\n";
    my $link_id = 1;
    foreach my $cM ( keys %{ $chr_cM_bounds{$chr} } ) {
        my $cM_start_pos = min( keys %{ $chr_cM_bounds{$chr}{$cM} } );
        my $cM_end_pos   = max( keys %{ $chr_cM_bounds{$chr}{$cM} } );
        my $scf1         = $chr_cM_bounds{$chr}{$cM}{$cM_start_pos};
        my $scf2         = $chr_cM_bounds{$chr}{$cM}{$cM_end_pos};
        next if ( $scf1 eq $scf2 );

        my $cM_physical_distance = $cM_end_pos - $cM_start_pos;
        print $linkage_xml "\t<CONTIG\n";
        print $linkage_xml "\t  ID=\"$scf1\"\n";
        print $linkage_xml "\t  NAME=\"$scf1\"\n";
        print $linkage_xml "\t>\n";
        print $linkage_xml "\t<CONTIG\n";
        print $linkage_xml "\t  ID=\"$scf2\"\n";
        print $linkage_xml "\t  NAME=\"$scf2\"\n";
        print $linkage_xml "\t>\n";

        print $linkage_xml "\t<LINK\n";
        print $linkage_xml "\t  ID=\"link_$link_id\"\n";
        print $linkage_xml "\t  SIZE=\"$cM_physical_distance\"\n";
        print $linkage_xml "\t  TYPE=\"Linkage\"\n";
        print $linkage_xml "\t>\n";
        $link_id++;
    }

    print $linkage_xml "</EVIDENCE>\n";
    close $linkage_xml;

    open my $info_file, '>', "$dir_name/$chr.info"
      or croak "Can't open info file for $chr $OS_ERROR\n";

    my %scf_reads;
    my %scf_bridges;
    foreach my $scf1 ( keys %{ $chr_scfs{$chr} } ) {
        print $info_file "Input\t$scf1\n";

        foreach my $ctg_id ( keys %{ $scf_contigs{$scf1} } ) {
            my $ctg_name = $scf_contigs{$scf1}{$ctg_id}{name};
            foreach my $read ( keys %{ $ctg_bridges{$ctg_name} } ) {

                next if ( defined $read_seen{$read} );

                $all_reads{$read}++;
                my $end1 = $ctg_bridges{$ctg_name}{$read};
                my $end2 = $end1 eq 'a' ? 'b' : 'a';

                if ( ( !defined $mates{$read}{ends}{$end2}{ctg} ) ) {
                    $missing_end_reads{$read}++;
                    next;
                }

                my $scf2 = $contigs_scf{ $mates{$read}{ends}{$end2}{ctg} }{scf};

                if ( $scf1 eq $scf2 ) {
                    $same_scf_reads{$read}++;
                    next;
                }

                if ( defined $scf_not_on_map{$scf2} ) {
                    $off_map_reads{$read}++;
                    next;
                }
                if ( !defined $chr_scfs{$chr}{$scf2} ) {
                    $off_chrom_reads{$read}++;
                    next;
                }

                $spanning_reads{$read}++;

                my $scf1_link_end =
                  $mates{$read}{ends}{$end1}{scf_orient} eq 'RC' ? 'B' : 'E';
                my $scf1_span =
                  ( $scf1_link_end eq 'B' )
                  ? $mates{$read}{ends}{$end1}{scf_start}
                  : $scf_lengths{$scf1} - $mates{$read}{ends}{$end1}{scf_end};

                my $scf2_link_end =
                  $mates{$read}{ends}{$end2}{scf_orient} eq 'RC' ? 'B' : 'E';
                my $scf2_span =
                  ( $scf2_link_end eq 'B' )
                  ? $mates{$read}{ends}{$end2}{scf_start}
                  : $scf_lengths{$scf2} - $mates{$read}{ends}{$end2}{scf_end};

                my $bridge_length = $scf1_span + $scf2_span;
                print $info_file
"Mate\t$chr\t$read\t$scf1\t$end1\t$scf_lengths{$scf1}\t$mates{$read}{ends}{$end1}{scf_start}\t$mates{$read}{ends}{$end1}{scf_end}\t$mates{$read}{ends}{$end1}{scf_orient}\t$scf1_span\t$scf2\t$end2\t$scf_lengths{$scf2}\t$mates{$read}{ends}{$end2}{scf_start}\t$mates{$read}{ends}{$end2}{scf_end}\t$mates{$read}{ends}{$end2}{scf_orient}\t$scf2_span\t$scf1-$scf1_link_end -- ";
                printf $info_file "%6d", $bridge_length;
                print $info_file " -- $scf2_link_end-$scf2";

                if ( ( $scf1_span > 8000 ) || ( $scf2_span > 8000 ) ) {
                    $distant_reads{$read}++;
                    print $info_file "\tDistant; rejected\n";
                    next;
                }

                print $info_file "\tAccepted\n";
                $scf_reads{$scf1}{$read} = $end1;
                $scf_reads{$scf2}{$read} = $end2;

                $scf_bridges{$scf1}{$scf2}{libraries}
                  { $mates{$read}{library} }++;
                $scf_bridges{$scf1}{$scf2}{reads}{$read}++;

                $read_seen{$read}++;
            }
        }

    }

    open my $fasta_file, '>', "$dir_name/$chr.fasta"
      or croak "Can't open FASTA output file for $chr $OS_ERROR\n";

    open my $contigs_file, '>', "$dir_name/$chr.contig"
      or croak "Can't open contigs output file for $chr $OS_ERROR\n";

    open my $mates_file, '>', "$dir_name/$chr.mates"
      or croak "Can't open mates output file for $chr $OS_ERROR\n";

    foreach my $library ( keys %library ) {
        print $mates_file
"library\t$library\t$library{$library}{min}\t$library{$library}{max}\n";

    }

    # Output mate pairs; filter by read count
    foreach my $scf1 ( keys %scf_bridges ) {
        foreach my $scf2 ( keys %{ $scf_bridges{$scf1} } ) {

            print $info_file "Link\t$chr\t$scf1\t$scf2";

            foreach
              my $library ( keys %{ $scf_bridges{$scf1}{$scf2}{libraries} } )
            {
                print $info_file
                  "\t$library=$scf_bridges{$scf1}{$scf2}{libraries}{$library}";
            }

            if ( scalar( keys %{ $scf_bridges{$scf1}{$scf2}{reads} } ) <= 1 ) {

                print $info_file "\tSingleton";

                # Check that the scaffolds share a cM locus
                my $shared_group = 0;
                foreach my $cM ( keys %{ $chr_scfs{$chr}{$scf1} } ) {
                    if ( defined $chr_scfs{$chr}{$scf2}{$cM} ) {
                        $shared_group++;
                        print $info_file "\t$cM";
                    }
                }
                if ( !$shared_group ) {
                    print $info_file "\tNo shared locus\tRejected\n";
                    next;
                }
                else {
                    print $info_file "\tShared locus";
                }
            }
            else {
                print $info_file "\tMultiple reads";
            }
            print $info_file "\tOK\n";

            foreach my $read ( keys %{ $scf_bridges{$scf1}{$scf2}{reads} } ) {
                my $end1  = $scf_reads{$scf1}{$read};
                my $end2  = $scf_reads{$scf2}{$read};
                my $mate1 = $end1 eq 'a' ? 1 : 2;
                my $mate2 = $end1 eq 'a' ? 2 : 1;
                print $mates_file
"$read$end1\.$mate1\t$read$end2\.$mate2\t$mates{$read}{library}\n";
            }
        }
    }

    my %readseq_checklist;
    foreach my $scf ( keys %scf_reads ) {

        # Write contig file header for this scaffold
        my $scf_read_num = scalar keys %{ $scf_reads{$scf} };
        print $contigs_file
          "##$scf $scf_read_num $scf_lengths{$scf} bases, 00000000 checksum.\n";

        # Construct scaffold sequence from contig sequences
        my $scf_seq = "";
        foreach my $ctg_id (
            sort { $a <=> $b }
            keys %{ $scf_contigs{$scf} }
          )
        {
            my $ctg_name = $scf_contigs{$scf}{$ctg_id}{name};
            my $ctg_seq =
              substr $contigs{ $scf_contigs{$scf}{$ctg_id}{name} },
              $contigs_scf{$ctg_name}{ctg_start} - 1,
              $contigs_scf{$ctg_name}{ctg_end};
            my $ctg_orient = $scf_contigs{$scf}{$ctg_id}{orient};
            if ( $ctg_orient eq '-' ) {
                $ctg_seq = reverse $ctg_seq;
            }
            $scf_seq .= $ctg_seq;
            if ( defined $scf_contigs{$scf}{$ctg_id}{gap} ) {
                for ( 1 .. $scf_contigs{$scf}{$ctg_id}{gap} ) {
                    $scf_seq .= "-";
                }
            }
        }

        # Write scaffold sequence to FASTA and contig files
        print $fasta_file ">$scf\n";
        print $fasta_file "$scf_seq\n";
        print $contigs_file "$scf_seq\n";

        foreach my $read ( keys %{ $scf_reads{$scf} } ) {
            my $read_ref    = $mates{$read}{ends}{ $scf_reads{$scf}{$read} };
            my $read_length = $read_ref->{scf_end} - $read_ref->{scf_start} + 1;
            my $read_seq    = substr $scf_seq, $read_ref->{scf_start},
              $read_length;
            my $mate = $scf_reads{$scf}{$read} eq 'a' ? 1 : 2;

            print $contigs_file "#$read$scf_reads{$scf}{$read}\.$mate("
              . ( $read_ref->{scf_start} - 1 )
              . ") [$read_ref->{scf_orient}] $read_length bases, 00000000 checksum. {";
            if ( $read_ref->{scf_orient} eq 'RC' ) {
                print $contigs_file "$read_length 1";
            }
            else { print $contigs_file "1 $read_length" }
            print $contigs_file
              "} <$read_ref->{scf_start} $read_ref->{scf_end}>\n$read_seq\n";
            if ( !defined $readseq_checklist{"$read$scf_reads{$scf}{$read}"} ) {
                print $fasta_file
                  ">$read$scf_reads{$scf}{$read}\.$mate\n$read_seq\n";
            }
        }
    }

    close $contigs_file;
    close $fasta_file;
    close $info_file;
    close $mates_file;

    my $output = `echo $chr > $dir_name/$chr.bambus.log`;
    print STDERR "\ttoAmos";
    $output =
`toAmos -c $dir_name/$chr.contig -s $dir_name/$chr.fasta -m $dir_name/$chr.mates -o $dir_name/$chr.afg 1>> $dir_name/$chr.bambus.log 2>&1`;
    print STDERR "\tbank-transact";
    $output =
`bank-transact -b $dir_name/$chr.bnk -m $dir_name/$chr.afg -cf 1>> $dir_name/$chr.bambus.log 2>&1`;
    print STDERR "\tclk";
    $output = `clk -b $dir_name/$chr.bnk 1>> $dir_name/$chr.bambus.log 2>&1`;
    print STDERR "\tBundler";
    $output =
      `Bundler -b $dir_name/$chr.bnk 1>> $dir_name/$chr.bambus.log 2>&1`;
    print STDERR "\tMarkRepeats";
    $output =
`MarkRepeats -b $dir_name/$chr.bnk > $dir_name/$chr.repeats 2>> $dir_name/$chr.bambus.log`;
    print STDERR "\tOrientContigs";
    $output =
`OrientContigs -b $dir_name/$chr.bnk -prefix $dir_name/$chr -repeats $dir_name/$chr.repeats -noreduce -redundancy 1 1>> $dir_name/$chr.bambus.log 2>&1`;
    print STDERR "\tDONE\n";
}

print STDERR "Mate pairs with one end on linkage map="
  . ( scalar keys %all_reads )
  . "\nMate pair ends mapping to contigs without scaffolds (haplotypes)="
  . ( scalar keys %missing_end_reads )
  . "\nMate pair ends mapping to scaffolds off linkage map="
  . ( scalar keys %off_map_reads )
  . "\nMate pair ends mapping to scaffolds on another chromosome="
  . ( scalar keys %off_chrom_reads )
  . "\nMate pairs mapping to same scaffold="
  . ( scalar keys %same_scf_reads )
  . "\nMate pairs spanning different scaffolds="
  . ( scalar keys %spanning_reads )
  . "\nMate pairs mapping >20 Kb apart="
  . ( scalar keys %distant_reads )
  . "\nUsable mate pairs="
  . ( scalar( keys %spanning_reads ) - scalar( keys %distant_reads ) ) . "\n";

sub find_link_end {
    my ( $orient, $end ) = @_;
    my $link_end;
    if ( $orient eq 'RC' ) {
        $link_end = ( $end eq 'a' ) ? 'E' : 'B';
    }
    else {
        $link_end = ( $end eq 'a' ) ? 'B' : 'E';
    }
    return $link_end;
}
