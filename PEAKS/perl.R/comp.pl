use strict; use warnings; 
my $longfl = 'long.txt';
my $shortfl = 'short.txt';

open(IN, "<$longfl" ); my @linesL = <IN>; close (IN);  
my %linesL = (); 
foreach my $line (@linesL) {
 chomp $line; 
 my ($chr, $pos, @res) = split( /\s+/, $line );
 $linesL{"chr-$pos"} = $line; 
 print "LONG[$chr $pos $line]\n"; 
}
#close(IN);

open(IN, "<$shortfl" ); 
my %linesS = (); 
while(my $line =<IN>) {
 chomp $line; 
 my ($chr, $pos, @res) = split( /\s+/, $line );
 $linesS{"chr-$pos"} = $line; 
 print "SHORT[$chr $pos $line]\n"; 
}
close(IN);

#find out unique short and long SNPs
open(OUT, ">_diffSNPS.txt"); 
foreach my $key (keys %linesS) {
  if (exists $linesL{$key}) { #do nothing
    #   print $linesL{$key}; 
  } else {
    print OUT $linesS{$key}." short\n"; 
 }
}

foreach my $key (keys %linesL) {
  if (exists $linesS{$key}) { #do nothing
  } else {
    print OUT $linesL{$key}." long\n"; 
 }
}

close(OUT);
