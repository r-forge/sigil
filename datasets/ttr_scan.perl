#!/usr/bin/perl
# -*-cperl-*-
## Scan CWB corpus to collect type/token counts and similar statistics for each text 
$| = 1;

use strict;
use warnings;

use CWB;
use CWB::CL::Strict;

die "Usage:  ttr_scan.perl CORPUS ttr_table.txt.gz [lemma]\n"
  unless @ARGV == 2 or @ARGV == 3;

our $Corpus = uc(shift @ARGV);
our $Outfile = shift @ARGV;
our $Attribute = (@ARGV) ? lc(shift @ARGV) : "word";

our $C = new CWB::CL::Corpus $Corpus;
our $W = $C->attribute($Attribute, "p");
our $POS = $C->attribute("pos", "p");
our $Text = $C->attribute("text_id", "s");

our %nominal = map {$_ => 1} $POS->regex2id("N.*"); # IDs of POS tags for nominals

our $OUT = CWB::OpenFile "> $Outfile";
print $OUT "tokens types hapaxes nouns wordlen avgfreq\n";

my $n_texts = $Text->max_struc;
my $old_perc = -1;
foreach my $t (0 .. $n_texts - 1) {
  my $id = $Text->struc2str($t);
  $id =~ s/\s+/_/g; $id =~ s/\\//g; $id =~ s/"/''/g; # remove dangerous characters from ID
  my ($start, $end) = $Text->struc2cpos($t);
  my @token_id = $W->cpos2id($start .. $end);
  my $n_tokens = @token_id;
  my $total_wordlen = sum($W->id2strlen(@token_id));
  my $total_freq = sum($W->id2freq(@token_id));
  my %types = ();
  foreach (@token_id) { $types{$_}++ }
  my $n_types = keys %types;
  my $n_hapax = grep { $_ == 1 } values %types;
  my $n_nominal = grep { $nominal{$_} } $POS->cpos2id($start .. $end);
  printf $OUT "\"%s\" %d %d %d %d %.6f %.3f\n", $id, $n_tokens, $n_types, $n_hapax, $n_nominal, $total_wordlen / $n_tokens, $total_freq / $n_tokens;
  my $perc = 100 * ($t+1) / $n_texts;
  if ($perc - $old_perc >= .01) {
    printf "%8d / %d  | %6.2f%s complete\r", $t+1, $n_texts, $perc, '%';
    $old_perc = $perc;
  }
}
print " " x 60, "\r";
$OUT->close;

print "Statistics for $n_texts texts written to file $Outfile\n";

# --------------------------------------------------

sub sum {
  my $sum = 0;
  foreach (@_) { $sum += $_ }
  return $sum;
}