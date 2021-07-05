#!/usr/bin/perl

my $line = <STDIN>;

print "$line";

my %hashLines = ();
while($line = <STDIN>){
    chomp($line);

    my @tokens = split(/,/,$line);

    $hashLines{$tokens[0]} = $line;
}

foreach $key (sort {$a <=> $b} keys %hashLines){
    print "$hashLines{$key}\n";
}
