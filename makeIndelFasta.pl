$usage = "<vcf> <fasta> [name] [indent] [lowfreq]";

open F, "<", shift or die $usage; 
$fastaFile = shift or die $usage;
$name = shift or $name = "";
$indent = shift or $indent = 25;
$lowFreq = shift or die $lowFreq = 0.25;

`mkdir -p tmp/pieces`;
open B, ">", "tmp/pieces/indel_pieces_$name.bed";
open O, ">", "tmp/pieces/indel_pieces_$name.inf";
open V, ">", "tmp/pieces/indel_pieces_$name.vcf"; # only selected places

print O "contig\tpos\taltFreq\tDP\tDP4\tseq\talt\n";

while(<F>){
	next if (/#/);
	@P = split(/\t/);
        
	# skip positions with more then 1 gaps
        next if(abs (length($P[4]) - length($P[3])) > 1);

	$start = $P[1] - $indent; # conserved letter
	next if($start < 0); 
	$end = $P[1] + $indent; 
	# TBD check contig size
	print B "$P[0]\t$start\t$end\t$P[0]____$P[1]\n";
	
	# info 
	my @inf = split(/:/,$P[9]);
	my $ref_dep = $inf[4];
	my $alt_dep = $inf[5];
	my $alt_freq = ($alt_dep+1)/($ref_dep+1);
	print O "$P[0]\t$P[1]\t$alt_freq\t$ref_dep\t$alt_dep\t$P[3]\t$P[4]\n"; 
	print V $_;
}

close (B);
close (O);
close (V);

