#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

&main;

sub main {
	my %opts = (t=>4, p=>'fmdef', k=>51, e=>29, l=>61);
	getopts('t:p:k:f:r:e:l:C', \%opts);
	die(qq/
Usage:   fm2-unitig.pl [options] <in.fq>\n
Options: -p STR     output prefix [$opts{p}]
         -k INT     min overlap length during unitig construction [$opts{k}]
         -l INT     min overlap length during graph cleaning [$opts{l}]
    	 -e INT     k-mer length used for error correction [$opts{e}]
         -t INT     number of threads [$opts{t}]
         -C         don't cut at low-quality bases
\n/) if (@ARGV == 0);

	$opts{f} ||= gwhich("fermi2");
	$opts{r} ||= gwhich("ropebwt2");
	die("[E::main] failed to find the 'fermi2' executable") unless (-x $opts{f});
	die("[E::main] failed to find the 'ropebwt2' executable") unless (-x $opts{r});

	my $is_file = (-f $ARGV[0]);

	my @lines = ();
	push(@lines, qq/PREFIX=$opts{p}/, '');
	push(@lines, qq/EXE_FERMI2=$opts{f}/, qq/EXE_ROPEBWT2=$opts{r}/);
	push(@lines, qq/K_UNITIG=$opts{k}/, qq/K_EC=$opts{e}/, qq/K_CLEAN=$opts{l}/);
	push(@lines, qq/N_THREADS=$opts{t}/, '');

	push(@lines, qq/all:\$(PREFIX).mag.gz/, "");

	push(@lines, qq/\$(PREFIX).raw.fmd:/);
	my $rb2_opts = defined($opts{C})? "-drq3" : "-drq20 -x31";
	if ($is_file) {
		push(@lines, qq/\t\$(EXE_ROPEBWT2) $rb2_opts $ARGV[0] > \$@ 2> \$@.log/);
	} else {
		push(@lines, qq/\t$ARGV[0] | \$(EXE_ROPEBWT2) $rb2_opts > \$@ 2> \$@.log/);
	}
	push(@lines, "");

	push(@lines, qq/\$(PREFIX).ec.fq.gz:\$(PREFIX).raw.fmd/);
	if ($is_file) {
		push(@lines, qq/\t\$(EXE_FERMI2) correct -t \$(N_THREADS) -k \$(K_EC) \$< $ARGV[0] 2> \$@.log | gzip -1 > \$@/);
	} else {
		push(@lines, qq/\t$ARGV[0] | \$(EXE_FERMI2) correct -t \$(N_THREADS) -k \$(K_EC) \$< \/dev\/stdin 2> \$@.log | gzip -1 > \$@/);
	}
	push(@lines, "");

	push(@lines, qq/\$(PREFIX).ec.fmd:\$(PREFIX).ec.fq.gz/);
	push(@lines, qq/\tgzip -dc \$< | \$(EXE_ROPEBWT2) -dNCr > \$@ 2> \$@.log/, "");

	push(@lines, qq/\$(PREFIX).flt.sub:\$(PREFIX).ec.fmd/);
	push(@lines, qq/\t\$(EXE_FERMI2) occflt -t \$(N_THREADS) \$< > \$@ 2> \$@.log/, "");

	push(@lines, qq/\$(PREFIX).flt.fmd:\$(PREFIX).ec.fmd \$(PREFIX).flt.sub/);
	push(@lines, qq/\t\$(EXE_FERMI2) sub -ct \$(N_THREADS) \$< \$(PREFIX).flt.sub > \$@ 2> \$@.log/, "");

	push(@lines, qq/\$(PREFIX).pre.mag.gz:\$(PREFIX).flt.fmd/);
	push(@lines, qq/\t\$(EXE_FERMI2) assemble -l \$(K_UNITIG) -t \$(N_THREADS) \$< 2> \$@.log | gzip -1 > \$@/, "");

	push(@lines, qq/\$(PREFIX).mag.gz:\$(PREFIX).pre.mag.gz/);
	push(@lines, qq/\t\$(EXE_FERMI2) simplify -CSo \$(K_CLEAN) \$< 2> \$@.log | gzip -1 > \$@/, "");

	print(join("\n", @lines), "\n");
}

sub gwhich {
    my $progname = shift;
    my $addtional_path = shift if (@_);
    my $dirname = &dirname($0);
    my $tmp;

    chomp($dirname);
    if ($progname =~ /^\// && (-x $progname)) {
        return $progname;
    } elsif (defined($addtional_path) && ($tmp = &which($progname, $addtional_path))) {
        return $tmp;
    } elsif (-x "./$progname") {
        return "./$progname";
    } elsif (defined($dirname) && (-x "$dirname/$progname")) {
        return "$dirname/$progname";
    } elsif (($tmp = &which($progname))) {
        return $tmp;
    } else {
        return;
    }
}

sub dirname {
	my $prog = shift;
	return '.' if (!($prog =~ /\//));
	$prog =~ s/\/[^\s\/]$//g;
	return $prog;
}
