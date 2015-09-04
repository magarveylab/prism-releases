#!/usr/bin/perl -l
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Species;
use Bio::Annotation::Collection;

my $seqin = Bio::SeqIO->new(-file   => $ARGV[0],
                                                        -format => 'genbank');

while (my $seq = $seqin->next_seq) {
        print $seq->desc ? $seq->desc: "";
        print $seq->accession ? $seq->accession: "";
        print join "\t", $seq->species->classification ?  $seq->species->classification : "";
        print $seq->species->genus ? $seq->species->genus: "";
        print $seq->species->species ? $seq->species->species: "";
        print $seq->species->sub_species ? $seq->species->sub_species: "" ;

        print "\n\n";
        for my $feat ($seq->get_SeqFeatures) {
                if ($feat->primary_tag eq "rRNA") {
                                #print $feat->primary_tag;
                        for my $tag ($feat->get_all_tags) {
                                if ($tag eq "product") {
                                        for my $value ($feat->get_tag_values($tag)) {
                                                if ($value =~ /16S ribosomal RNA/i) {
                                                    print "16S";
                                                        my $location = $feat->location;
                                                        print $location->start;
                                                        print $location->end;
                                                        print $location->strand;
                                                        for my $tag ($feat->get_all_tags) {
                                                                print $tag, "\t", $feat->get_tag_values($tag);
                                                                #if ($tag eq "translation") {
                                                                #       print "16S\t", $feat->get_tag_values($tag); 
                                                                #}
                                                        }
                                                        print "###";
                                                }
                                        }
                                }
                        }
                }
        }
}
