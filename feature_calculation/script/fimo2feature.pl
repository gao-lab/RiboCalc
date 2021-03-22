use warnings;
use strict;

open LIST,"<$ARGV[0]" or die $!; #fasta
open IN,"<$ARGV[1]" or die $!; #fimo output

my @motif = ("pwm_atg", "pwm_ctg", "pwm_gtg", "pwm_ttg");

# record the transcript list
my (@arr,$id,%hash);
while(<LIST>){
	chomp;
	if ($_=~/^>/){
		$id = (split,$_)[0];
		$id =~ s/^>//;
		for (@motif){
			$hash{$id}{$_} = 0;
		}
	}
}
close LIST;

# summary the fimo output
while(<IN>){
	chomp;
	next if ($_ =~ /^#/);
	@arr = split;
	next if ($arr[4] eq "-"); #only count for the sense strand
	if (exists $hash{$arr[1]}){
		$hash{$arr[1]}{$arr[0]} += 1 if (exists $hash{$arr[1]}{$arr[0]});
	}else{
		next;
	}
}
close IN;

my $head = join("\t",@motif);
print "ID\t$head\n";
for (keys %hash){
	$id = $_;
	print "$id\t";
	for (@motif){
		$hash{$id}{$_} = 1 if ($hash{$id}{$_} > 0); #use binary code for describing motif existance
		print "\t$hash{$id}{$_}";
	}
	print "\n";
}
