#!/usr/bin/perl

use Bio::Search::Result::BlastResult;
use Bio::SearchIO;
use List::Util qw[min max];
use Bio::SeqIO;
use v5.10;

$usage = "[blast xml] [vcf] [assembly fasta] [out name] [indent]";
$blastFile = shift or die $usage;
$vcfFile = shift or die $usage;
$fastaFName = shift or die $usage;
$outPrefix = shift or die $usage;
$indent = shift or die $usage;

# l o a d    v c f
open V, "<", $vcfFile; 
while(<V>){
	@P = split(/\t/);
	$VCF{$P[0]}{$P[1]}{"seq"} = $P[3];
	$VCF{$P[0]}{$P[1]}{"alt"} = $P[4];	
	$VCF{$P[0]}{$P[1]}{"contig"} = $P[0];
	$VCF{$P[0]}{$P[1]}{"pos"} = $P[1];	
}
close V;

# l o a d    a s s e m b l y   f a s t a 
my $assemIO = Bio::SeqIO->new( -file => "<$fastaFName", -format => "fasta");
while(my $seqobj = $assemIO->next_seq) {
    my $id  = $seqobj->display_id;    
    my $seq = uc $seqobj->seq;          
    $assembly{$id} = $seq;
}

# M a i n 
my %toCorrect;
my %indent; 
my $report = Bio::SearchIO->new( -file=>$blastFile, -format => blastxml);
`mkdir -p log`;
open L, ">log/$outPrefix.log";

while(my $result = $report->next_result){
	$qlen = $result->query_length; 
	$qname = $result->query_description; 
	@qname = split(/____/, $qname);
	my $contig = $qname[0];
	my $pos = $qname[1];
	
	# g e t    b l a s t   r e s u l t
	#  delete variants with no hits
	## look only first hit and hsp
	my $hit = $result->next_hit; 
	next if(!$hit);
	my $hsp = $hit->next_hsp; 
	next if(!$hsp);
	$hlen = $hsp->length('query');
	
	my $qstart = $hsp->start('query');
	my $qend = $hsp->end('query');

	# skip if alignment is too short
	if(abs($qend - $qstart) < 2*$indent -1) { print "qend-qstart:".($qend - $astart)."(int: $intend)" . "\n"; next; }; # CHECK!!!! 
	# info
	$gaps = $hsp->gaps;
	$evalue = $hsp->evalue;	
	$frac_identical = $hsp->frac_identical;
	
	# c o m p a r e    w i t h    v c f
	my $start = $indent-$qstart; 
	my $toend = $hit->length('total') - $start;
	my $seq = uc $VCF{$qname[0]}{$qname[1]}{"seq"};  
	my $alt =  uc $VCF{$qname[0]}{$qname[1]}{"alt"}; 
	$alt =~ s/,.*//; 	# select only first alt variant
	my $ref = uc substr($hsp->hit_string, $start, $toend);
	# print
	say L "- - -\n$qname"; 
 	say L "assembly:  " . substr ($hsp->query_string, $start, $toend);
	say L "blast hit: " . $ref;
	say L "seq_var:   " . $seq;
	say L "alt_var:   " . $alt;
	#say L "fasta_seq: " . substr($assembly{$contig}, $pos - 1, length($seq));
	# compare and remember
	$ref =~ s/-//g;
	$alt = uc $alt;
	if(substr($ref, 0, length($alt)) eq $alt) {
		next if(length($alt) < length($seq) and substr($ref, 0, length($seq)) eq $seq); 
		$dLength = length($alt) - length($seq);
		my $pos_original = $pos;
		$pos += $indent{$contig} - 1;
		next if($pos < 0);
		print L "pos: $pos\tindent: $indent{$contig}\n"; 
		$toCorrect{$contig}{$pos}{"seq"} = $seq; 
		$toCorrect{$contig}{$pos}{"alt"} = $alt; 
		$toCorrect{$contig}{$pos}{"evalue"} = $evalue;
		$toCorrect{$contig}{$pos}{"frac_ident"} = $frac_identical;
	
	$toCorrect{$contig}{$pos}{"pos"} = $pos_original;
		$indent{$contig} += $dLength;
		print L "will be corrected...\n";	
	}
}


# c o r r e c t

print L "CORRECTING...\n";
open C, ">log/${outPrefix}_corrected.txt";

my $corrections;
foreach my $contig (sort keys %toCorrect){
	foreach my $pos (sort { $a <=> $b } keys %{ $toCorrect{$contig} }) {
		# $pos indenting because of del/inserts 
		$seq = $toCorrect{$contig}{$pos}{"seq"};
		$alt = $toCorrect{$contig}{$pos}{"alt"};
		$evalue = $toCorrect{$contig}{$pos}{"evalue"};
		$frac_ident =$toCorrect{$contig}{$pos}{"frac_ident"};  
		$pos_original =$toCorrect{$contig}{$pos}{"pos"}; 
		$dLength = length($alt) - length($seq);
		substr($assembly{$contig}, $pos, length($seq)) = $alt;
		say L "seq: " . substr($assembly{$contig}, $pos, length($seq));
		say L "alt: $alt";
		say L "pos: $pos";
		say L "con: $contig";
		say L "---";
		say C "$contig\t$pos_original\t$evalue\t$frac_ident";
		$corrections++;
	}
}

close C;
# w r i t e 

my $fastaString;
foreach my $contig (sort keys %assembly){
	$fastaString .= ">$contig\n$assembly{$contig}\n";
}	

open($stringfh, "<", \$fastaString) or die "Could not open string for reading: $!";  
my $inseqIO = Bio::SeqIO->new(-fh => $stringfh, -format => "fasta");
my $outseqIO = Bio::SeqIO->new(-file => ">corrected/$outPrefix.fasta", -format => "fasta");
while (my $inseq = $inseqIO->next_seq) {
       $outseqIO->write_seq($inseq);
       }

print L "\n\n\n corrected: $corrections";

close(L);





