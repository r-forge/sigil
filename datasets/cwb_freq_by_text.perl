#!/usr/bin/perl
## -*-cperl-*-
## Run CQP queries on specified corpora, generating R table with frequency information per document
##
$| = 1;

use strict;
use warnings;

use CWB;
use CWB::CQP;
use Time::HiRes qw(time);
use Getopt::Long;

our $USAGE = <<'====';
Usage:  brown_lob_query.perl [options] tbl/infinitives.tbl split.inf='"to"%c [pos="RB.*"] [pos="V."]' ...
Options:
  --corpora=<list>,     comma-separated list of corpora to be scanned, use "label:ID" to prepend text IDs with "label_"_
      -D <list>             [default: brown:TT-BROWN,lob:TT-LOB]
  --id=<att>, -i <att>  s-attribute specifying text IDs [default: text_id]
  --sentences, -s       sentence counts (i.e. number of sentences containing matches)
====

our $Corpora = "brown:TT-BROWN,lob:TT-LOB"; # use Marco's TreeTagged version of BROWN and LOB -> same attributes and tagsets
our $ID_ATT = "text_id";
our $SentenceCounts = 0;
our $Help = 0;

my $ok = GetOptions(
  "corpora|D=s" => \$Corpora,
  "id|i=s" => \$ID_ATT,
  "sentences|sent|s" => \$SentenceCounts,
  "help|h" => \$Help,
);

die $USAGE
  unless @ARGV >= 2 and $ok and not $Help;

our $OutFile = shift @ARGV;

our @Queries = ();
foreach my $item (@ARGV) {
  $item =~ /^([A-Za-z0-9._-]+)=(.+)$/
    or die "Syntax error in ``$item'' -- expected <name>='<query>'\n";
  push @Queries, {-name => $1, -query => $2};
}

our @Corpora = map {
  die "Format error in corpus specified ``$_'' -- expected [<label>:]<ID>\n"
    unless /^(?:([A-Za-z0-9_-]+):)?([A-Z0-9_-]+)$/;
  (defined $1) ? ["$1_", $2] : ["", $2];
} split /\s*,\s*/, $Corpora;

our $FH = CWB::OpenFile "> $OutFile";
print $FH join("\t", "id", map {$_->{-name}} @Queries), "\n";

foreach my $iteration (@Corpora) {
  our ($Label, $CORPUS) = @$iteration;

  our $CQP = new CWB::CQP;
  $CQP->set_error_handler('die');
  $CQP->exec($CORPUS); 

  print "Getting text IDs for $CORPUS ... ";
  $CQP->exec("IDs = <$ID_ATT> []");
  our @Texts = $CQP->exec("tabulate IDs match $ID_ATT");
  our $n_texts = @Texts;
  print "$n_texts texts\n";

  our %F = ();
  foreach my $Q (@Queries) {
    my $name = $Q->{-name};
    my $query = $Q->{-query};
    print "Searching '$name' ... ";
    my $T0 = time;
    my $expand = ($SentenceCounts) ? "within s expand to s" : "";
    $CQP->exec("A = $query $expand");
    my ($n) = $CQP->exec("size A");
    print "tabulating ... ";
    $CQP->run("tabulate A match $ID_ATT");
    while (my @row = $CQP->getrow) {
      $F{$name}{$row[0]}++;  # increment frequency count for text ID
    }
    my $dT = time - $T0;
    printf "%d matches in %.1f s\n", $n, $dT;
  }
  undef $CQP;

  my $T0 = time;
  print "Storing results ... ";
  foreach my $id (@Texts) {
    print $FH $Label.$id;
    foreach my $name (map {$_->{-name}} @Queries) {
      my $f = $F{$name}{$id} || 0;
      print $FH "\t", $f;
    }
    print $FH "\n";
  }
  my $dT = time - $T0;
  printf "%d lines in %.1f s\n", $n_texts, $dT;
}

print "Results have been saved to '$OutFile'.\n";
$FH->close;











