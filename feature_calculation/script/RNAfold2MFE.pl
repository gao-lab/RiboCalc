use warnings;
use strict;

open IN,"<$ARGV[0]" or die $!;
open OUT,">$ARGV[1]" or die $!;

my ($id,$line,$struct,$mfe);
print OUT "ID\tinit_fold\n";
while(<IN>){
	chomp;
	if ($_=~/^>/){
		$_ =~ s/^>//;
		$id = (split,$_)[0];
		<IN>;
		$line = <IN>;
		$mfe = (split/\s+/,$line)[-1];
		$mfe =~ s/^\(//;
		$mfe =~ s/\)$//;
		$mfe = abs($mfe);
		print OUT "$id\t$mfe\n";
	}
}
close IN;
close OUT;
