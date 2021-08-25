#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

open IN,"<$ARGV[0]" || die $!;
my $fasta = $ARGV[1];
my $ina = Bio::SeqIO->new(-file => $fasta, -format => 'fasta');
my (@arr,%hash,$id,$desc,$seq,$tid,@ids,$tpm);
while (my $obj = $ina->next_seq()){
        $id = $obj->id;
        $desc = $obj->desc;
        $seq = $obj->seq;
	$tid = (split /\|/, $id)[0];
	chomp($seq);
	unless ($seq eq "Sequenceunavailable" or $seq eq ""){
		$hash{$tid} = $seq;
	}
}

open OUT,">$ARGV[2]" || die $!;
open FA,">$ARGV[3]" || die $!;
while (<IN>){
	chomp;
	@arr = split;
	@ids = split /\|/, $arr[0];
	$tid = (split /\./, $ids[0])[0];
	$tpm = $arr[81]; #the column number for RiboTPM
	if (exists $hash{$tid}){
		print FA ">$tid\n$hash{$tid}\n";
		print OUT "$tid\t$tpm\n";
	}
}
close IN;
close OUT;
close FA;
