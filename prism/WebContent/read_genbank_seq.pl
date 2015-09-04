#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Species;
use Bio::Annotation::Collection;

my $seqin = Bio::SeqIO->new(-file   => $ARGV[0],
                                                      -format => 'genbank');
while (my $seq = $seqin->next_seq) {
        print $seq->seq ? $seq->seq : "NOSEQ";
}
