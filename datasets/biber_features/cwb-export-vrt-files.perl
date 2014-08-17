#!/usr/bin/perl
# -*-cperl-*-
## Export CWB-encoded corpus into separate .vrt file for each text
$| = 1;

use strict;
use warnings;

use CWB;
use CWB::CL::Strict;

die "Usage:  cwb-export-vrt-files.perl <CORPUS> <prefix> <id_att> <chunk_att> <att1> <att2> ...\n"
  unless @ARGV >= 5;

our ($Corpus, $Prefix, $IDAtt, $ChunkAtt, @Attributes) = @ARGV;

our $C = new CWB::CL::Corpus $Corpus;
our $TextID = $C->attribute($IDAtt, "s");
our $Chunk = $C->attribute($ChunkAtt, "s");
our @A = map { $C->attribute($_, "p") } @Attributes;

my $n_texts = $TextID->max_struc;
my $n_tokens = 0;

foreach my $n (0 .. $n_texts - 1) {
  my $id = $TextID->struc2str($n);
  my ($start, $end) = $TextID->struc2cpos($n);

  my ($first_chunk, $last_chunk) = $Chunk->cpos2struc($start, $end);
  die "Chunk <$ChunkAtt> doesn't align with start of text $id.\n"
    unless ($Chunk->struc2cpos($first_chunk))[0] == $start;
  die "Chunk <$ChunkAtt> doesn't align with end of text $id.\n"
    unless ($Chunk->struc2cpos($last_chunk))[1] == $end;

  my $OutFile = $Prefix . $id . ".txt";
  my $FH = CWB::OpenFile ">", $OutFile;  

  foreach my $c ($first_chunk .. $last_chunk) {
    printf $FH "<%s>\n", $ChunkAtt;
    ($start, $end) = $Chunk->struc2cpos($c);
    my $len = $end - $start + 1;
    my @cpos = $start .. $end;
    my @A_tokens = map { [$_->cpos2str(@cpos)] } @A;
    foreach my $i (0 .. $len - 1) {
      my @line = ();
      foreach my $att (@A_tokens) {
        push @line, $att->[$i];
      } 
      print $FH join("\t", @line), "\n";
    }
    printf $FH "</%s>\n", $ChunkAtt;
    $n_tokens += $len;
  }
  
  close $FH;
  printf "%8d / %d texts exported  %6.1fM tokens  \r", $n + 1, $n_texts, $n_tokens / 1e6
    if ($n & 0xf) == 0;
}

print " " x 60, "\r";
printf "%d texts exported to %s*.vrt, with %.1fM tokens\n", $n_texts, $Prefix, $n_tokens / 1e6;

