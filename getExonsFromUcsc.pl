#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use DBI;
use LWP::Simple;
use Pod::Usage;
use Data::Dumper;
use FindBin qw($RealBin);
use lib "$RealBin/lib/dapPerlGenomicLib";
use SortCoordinates;

my @genes;
my %opts = 
(
    genes => \@genes,
    build => "hg19",
);

GetOptions
(
    \%opts,
    "genes|g=s{,}",
    'l|list=s',
    "build|b=s",
    "flanks|f=i",
    "coding|c",
    "merge|m",
    "keep_info|k",
    "w|whole_gene",
    "s|standard_chr_only",
    "failures=s",
    "help|h",
    "manual",
) or pod2usage(-exitval => 2, -message => "Syntax error");
pod2usage (-verbose => 2) if $opts{manual};
pod2usage (-verbose => 1) if $opts{help};
pod2usage
(
    -exitval => 2, 
    -message => "--genes argument required - for help use --help"
) if not $opts{l} and not @genes;

if ($opts{l}){
    open (my $GENES, $opts{l}) or die "Could not open --gene_list '$opts{l}' for reading: $!\n";
    while (my $line = <$GENES>){
        my @s = split(/\s+/, $line); 
        push @genes, $s[0];
    }
}
if (not @genes){
    die "No gene IDs identified!\n";
}
$opts{flanks} = 0 if not defined $opts{flanks}; #default flanks value
my $FAILS_FH = undef;
if ($opts{failures}){
    open ($FAILS_FH, ">", $opts{failures}) or die "Could not create ".
            "$opts{failures} file: $!\n";
}

#CONNECT AND EXECUTE QUERY
my $dbh = DBI->connect_cached
(
    "dbi:mysql:$opts{build}:genome-mysql.cse.ucsc.edu", 
    "genome", 
    ''
);
my @regions = ();
if ($opts{w}){
    @regions = get_gene()
}else{
    @regions = get_exons();
}

#SORT AND MERGE OUR EXONS+FLANKS ARRAY
@regions = SortCoordinates::sortByCoordinate(array => \@regions);
if ($opts{merge}){
    @regions = SortCoordinates::mergeByCoordinate
    (
        array     => \@regions,
        keep_info => $opts{keep_info},
    );
}

foreach my $exon (@regions){
    if ($opts{s}){ #skip haplotype chromosomes
        next if $exon =~ /^(chr)?[0-9XYMT]+_/ or $exon =~ /^(chr)?Un_/;
    }
    print "$exon\n";   
}

if ($FAILS_FH){
    close $FAILS_FH;
}

#################################################
sub get_gene{
    my @tx_regions = ();
    my @cds_regions = ();
    my $command = 
    "SELECT name, chrom, txStart, txEnd, cdsStart, cdsEnd, strand, name2"
    ." FROM $opts{build}.refGene WHERE name2=?";
    my $sth = $dbh->prepare($command);
    foreach my $gene (@genes){
        my $found = 0;
        $sth->execute($gene);
        #PARSE RESULTS INTO EXON + FLANK REGIONS AND CDS REGIONS
        while (my  @row = $sth->fetchrow_array ) {
            $found++;
            push @cds_regions, join
            (
                "\t", 
                (
                    $row[1], 
                    $row[4] - $opts{flanks}, 
                    $row[5] + $opts{flanks},
                    "$gene|$row[0]",
                    0,
                    $row[6],
                )
            );
            push @tx_regions, join
            (
                "\t", 
                (
                    $row[1], 
                    $row[2] - $opts{flanks}, 
                    $row[3] + $opts{flanks},
                    "$gene|$row[0]",
                    0,
                    $row[6],
                )
            );
            
        }
        if (not $found){
            warn "No results for '$gene'\n";
            if ($FAILS_FH){
                print $FAILS_FH "$gene\n";
            }
        }
    }
    if ($opts{coding}){
        return @cds_regions;
    }else{
        return @tx_regions;
    }
}

#################################################
sub get_exons{
    my @cds_exons = ();#bed line for each cds region
    my @all_exons = ();#bed line for each exon and flanks
    my $command = 
    "SELECT name, chrom, cdsStart, cdsEnd, exonStarts, exonEnds, strand, name2"
    ." FROM $opts{build}.refGene WHERE name2=?";  
    my $sth = $dbh->prepare($command);
    foreach my $gene (@genes){
        my $found = 0;
        $sth->execute($gene);
        #PARSE RESULTS INTO EXON + FLANK REGIONS AND CDS REGIONS
        while (my  @row = $sth->fetchrow_array ) {
            $found++;
            my @starts = split(",", $row[4]);
            my @ends = split(",", $row[5]);
            #foreach exon check if coding and if so put in @cds_exons
            #add a region of exon + $opts{flanks} to each end to @all_exons
            for (my $i = 0; $i < @starts; $i++){
                my $ex;
                if ($row[6] eq '-'){
                    $ex = scalar(@starts) - $i ;
                }else{
                    $ex = $i + 1;
                }
                if ($starts[$i] >= $row[2] && $starts[$i] <= $row[3]){#after cds start and before cds end
                    if ($ends[$i] < $row[3]){#end is before cds end
                        push @cds_exons, join
                        (
                            "\t", 
                            (
                                $row[1], 
                                $starts[$i]- $opts{flanks}, 
                                $ends[$i]+ $opts{flanks},
                                "$gene|$row[0]_ex$ex",
                                0,
                                $row[6],
                            )
                                 
                        );
                    }else{#cds end in same exon
                        push @cds_exons, join
                        (
                            "\t", 
                            (
                                $row[1], 
                                $starts[$i]- $opts{flanks}, 
                                $row[3]+ $opts{flanks},
                                "$gene|$row[0]_ex$ex",
                                0,
                                $row[6],
                            )
                        );
                    }
                }elsif($starts[$i] < $row[2] and $row[2] < $ends[$i]){#cds start in middle of exon
                    if ($ends[$i] < $row[3]){#end is before cds end
                        push @cds_exons, join
                        (
                            "\t", 
                            (
                                $row[1], 
                                $row[2]- $opts{flanks}, 
                                $ends[$i]+ $opts{flanks},
                                "$gene|$row[0]_ex$ex",
                                0,
                                $row[6],
                            )
                        );
                    }else{#cds end in same exon
                        push @cds_exons, join
                        (
                            "\t", 
                            (
                                $row[1], 
                                $row[2]- $opts{flanks}, 
                                $row[3]+ $opts{flanks},
                                "$gene|$row[0]_ex$ex",
                                0,
                                $row[6],
                            )
                        );
                    }
                }
                    

                push @all_exons, join
                (
                    "\t", 
                    (
                        $row[1], 
                        $starts[$i] - $opts{flanks}, 
                        $ends[$i] + $opts{flanks},
                        "$gene|$row[0]_ex$ex",
                        0,
                        $row[6],
                    )
                );
            }
        }
        if (not $found){
            warn "No results for '$gene'\n";
            if ($FAILS_FH){
                print $FAILS_FH "$gene\n";
            }
        }
    }
    if ($opts{coding}){
        return @cds_exons;
    }else{
        return @all_exons;
    }
}

#################################################

=head1 NAME

getExonsFromUcsc.pl - get RefSeq exons for a given gene symbol from UCSC

=head1 SYNOPSIS

        getExonsFromUcsc.pl -g [gene symbol]  [options]
        getExonsFromUcsc.pl -h (display help message)
        getExonsFromUcsc.pl -m (display manual page)


=cut

=head1 ARGUMENTS

=over 8

=item B<-g    --genes>

Gene symbol to search for.

=item B<-b    --build>

Genome build/version to use. Default = hg19.

=item B<-f    --flanks>

Amount of flank bp to add to each exon (default = 0).

=item B<-m    --merge>

Merge overlapping exons into single regions.

=item B<-k    --keep_info>

When merging overlapping exons, use this option to print a column of details 
from each merged region.

=item B<-c    --coding>

Only output coding regions. 

=item B<-w    --whole_gene>

Output transcription start and stop coordinates only for each gene, or if 
--coding option is used, the CDS start and end coordinates.

=item B<-s    --standard_chr_only>

Do not output regions on haplotype or unplaced contigs.

=item B<-h    --help>

Display help message.

=item B<--manual>

Show manual page.


=back 

=cut


=head1 DESCRIPTION

This program outputs exon coordinates in BED format for given genes. The genome
version defaults to hg19, but can be changed via the -b/--build option.

=head1 EXAMPLES

        getExonsFromUcsc.pl -g ABCD1 > ABCD1_exons.bed
        
        getExonsFromUcsc.pl -g ABCD1 -f 100 > ABCD1_exons_plus100.bed

        getExonsFromUcsc.pl -g ABCD1 -b mm9 > mouse_ABCD1_exons.bed
        
        getExonsFromUcsc.pl -g ABCD1 ABCD2 > ABCD1_and_ABCD2_exons.bed

        getExonsFromUcsc.pl -g COL13A1 -c > COL13A1_coding_exons.bed
        
        getExonsFromUcsc.pl -g COL13A1 -c -m > COL13A1_coding_exons_merged.bed
        
        getExonsFromUcsc.pl -g COL13A1 -c -m -k > COL13A1_coding_exons_merged_with_info.bed



=cut

=head1 AUTHOR

David A. Parry

University of Edinburgh

=head1 COPYRIGHT AND LICENSE

Copyright 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut


