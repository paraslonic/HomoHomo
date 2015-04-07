open F, "<", shift;

$name_prev = "";

$str_prev = "";


while(<F>){
	next if(!/>/);
	@P = split(/\|/);
	$name = $P[2];
	if($name eq $name_prev and $name !~ /hypothetical/){
		print;	
	}	
	$name_prev = $name;
	$str_prev = $_;
}
