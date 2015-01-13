#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

die (qq/
Usage:   fermi2.pl <command> [arguments]\n
Command: unitig     generate Makefile for unitig assembly
         utglog     analyze log files generated by unitig
         mag2fmr    create FMR for multiple MAG unitig assemblies
\n/) if @ARGV == 0;

my $cmd = shift(@ARGV);
if ($cmd eq 'unitig') { &unitig(); }
elsif ($cmd eq 'mag2fmr') { &mag2fmr(); }
elsif ($cmd eq 'utglog') { &utglog(); }
else { die("ERROR: unknown command\n"); }

sub mag2fmr {
	my %opts = (l=>102, d=>3);
	getopts('ai:s:r:l:d:', \%opts);
	die (qq/fermi2.pl mag2fmr [-a] [-i in.fmr] <file1.mag.gz> [...]\n/) if @ARGV == 0;

	$opts{s} ||= gwhich("seqtk");
	$opts{r} ||= gwhich("ropebwt2");
	die ("ERROR: failed to find seqtk and ropebwt2") unless (-x $opts{s} && -x $opts{r});
	my @lines = ();
	my $prev = defined($opts{i})? $opts{i} : '';
	for my $fn (@ARGV) {
		unless (-f $fn) {
			warn("WARNING: skip non-existing file '$fn'");
			next;
		}
		my $pre = $fn =~ /(\S+)\.mag\.gz$/? $1 : $fn;
		push(@lines, qq/$pre.fmr:$fn $prev/);
		my $opt_rb2 = $prev? "-bRLi $prev" : "-bRL";
		my $tmp = qq/awk 'NR%2==0'|rev|tr "ACGT" "TGCA"|sort -S15G|tr "ACGT" "TGCA"|rev|$opts{r} $opt_rb2 > \$@ 2> \$@.log/;
		if (!defined($opts{a})) {
			my $genfa = qq/$opts{s} seq -nn -aq$opts{d} -l60 \$< | $opts{s} cutN -n1 - | $opts{s} seq -L$opts{l} -l0 | gzip -1 > $pre.fa.gz/;
			my $seqs = qq/(gzip -dc $pre.fa.gz; $opts{s} seq -rl0 $pre.fa.gz)/;
			push(@lines, qq/\t$genfa; $seqs|$tmp; rm -f $pre.fa.gz/, "");
		} else {
			my $seqs = qq/($opts{s} seq -l0 \$<; $opts{s} seq -rl0 \$<)/;
			push(@lines, qq/\t$seqs|$tmp/, "");
		}
		$prev = "$pre.fmr";
	}
	unshift(@lines, "all:$prev\n");

	print(join("\n", @lines), "\n");
}

sub unitig {
	my %opts = (t=>4, p=>'fmdef', l=>101, k=>-1, e=>0, T=>51, o=>-1, m=>-1, s=>'100m');
	getopts('t:p:k:f:r:c:e:l:m:CKEbs:T:', \%opts);
	die (qq/
Usage:   fermi2.pl unitig [options] <in.fq>\n
Options: -p STR     output prefix [$opts{p}]
         -s STR     approximate genome size (for bfc) [$opts{s}]
         -l INT     primary read length [$opts{l}]
         -T INT     use INT-mer for post-trimming\/filtering [$opts{T}]
         -k INT     min overlap length during unitig construction [based on -l]
         -o INT     min overlap length during graph cleaning [based on -l]
         -m INT     min overlap length for unambiguous merging [based on -l]
         -t INT     number of threads [$opts{t}]

         -e INT     fermi2 ec k-mer length (0 to use bfc) [$opts{e}]
         -K         don't drop reads during error correction
         -C         don't cut at low-quality bases for raw reads
         -E         cut at low-quality bases for corrected reads
\n/) if (@ARGV == 0);

	my $use_bfc = ($opts{e} == 0);

	if (!$use_bfc) {
		delete($opts{K}); delete($opts{C}); delete($opts{E});
	}

	$opts{k} = int($opts{l} * .5)  + 1 if $opts{k} < 0;
	$opts{m} = int($opts{l} * .75) + 1 if $opts{m} < 0;
	$opts{o} = $opts{k} + 5 if $opts{o} < 0;

	$opts{f} ||= gwhich("fermi2");
	$opts{r} ||= gwhich("ropebwt2");
	$opts{c} ||= gwhich("bfc");
	die("[E::main] failed to find the 'fermi2' executable") unless (-x $opts{f});
	die("[E::main] failed to find the 'ropebwt2' executable") unless (-x $opts{r});
	die("[E::main] failed to find the 'bfc' executable") unless (-x $opts{c});

	my $is_file = (-f $ARGV[0]);

	my @lines = ();
	push(@lines, qq/PREFIX=$opts{p}/, '');
	push(@lines, qq/EXE_FERMI2=$opts{f}/, qq/EXE_ROPEBWT2=$opts{r}/, qq/EXE_BFC=$opts{c}/);
	push(@lines, qq/K_UNITIG=$opts{k}/, qq/K_EC=$opts{e}/, qq/K_CLEAN=$opts{o}/, qq/K_TRIM=$opts{T}/, qq/K_MERGE=$opts{m}/, qq/GENOME_SIZE=$opts{s}/);
	push(@lines, qq/N_THREADS=$opts{t}/, "");
	push(@lines, qq/INPUT=$ARGV[0]/, "") unless ($is_file);

	push(@lines, qq/all:\$(PREFIX).mag.gz/, "");

	if ($use_bfc) {
		push(@lines, qq/\$(PREFIX).ec.fq.gz:/);
		if ($is_file) {
			push(@lines, qq/\t\$(EXE_BFC) -s \$(GENOME_SIZE) -t \$(N_THREADS) $ARGV[0] 2> \$@.log | gzip -1 > \$@/);
		} else {
			push(@lines, qq/\tbash -c '\$(EXE_BFC) -s \$(GENOME_SIZE) -t \$(N_THREADS) <(\$(INPUT)) <(\$(INPUT)) 2> \$@.log | gzip -1 > \$\@'/);
		}
		push(@lines, "");

		push(@lines, qq/\$(PREFIX).flt.fa.gz:\$(PREFIX).ec.fq.gz/);
		push(@lines, qq/\t\$(EXE_BFC) -1Qs \$(GENOME_SIZE) -k \$(K_TRIM) -t \$(N_THREADS) \$< 2> \$@.log | gzip -1 > \$@/, "");

		push(@lines, qq/\$(PREFIX).flt.fmd:\$(PREFIX).flt.fa.gz/);
		push(@lines, qq/\t\$(EXE_ROPEBWT2) -dNCr \$< > \$@ 2> \$@.log/, "");
	} else {
		push(@lines, qq/\$(PREFIX).raw.fmd:/);
		my $opt_rb2 = defined($opts{C})? "-drq3" : '-drq20 -x `expr $(K_EC) + 2`';
		if ($is_file) {
			push(@lines, qq/\t\$(EXE_ROPEBWT2) $opt_rb2 $ARGV[0] > \$@ 2> \$@.log/);
		} else {
			push(@lines, qq/\t$ARGV[0] | \$(EXE_ROPEBWT2) $opt_rb2 > \$@ 2> \$@.log/);
		}
		push(@lines, "");

		my $opt_ec = defined($opts{K})? "" : "-D";
		push(@lines, qq/\$(PREFIX).ec.fq.gz:\$(PREFIX).raw.fmd/);
		if ($is_file) {
			push(@lines, qq/\t\$(EXE_FERMI2) correct $opt_ec -t \$(N_THREADS) -k \$(K_EC) \$< $ARGV[0] 2> \$@.log | gzip -1 > \$@/);
		} else {
			push(@lines, qq/\t$ARGV[0] | \$(EXE_FERMI2) correct $opt_ec -t \$(N_THREADS) -k \$(K_EC) \$< \/dev\/stdin 2> \$@.log | gzip -1 > \$@/);
		}
		push(@lines, "");

		$opt_rb2 = defined($opts{E})? '-dCrq20 -x `expr $(K_UNITIG) + 2`' : '-dNCr';
		push(@lines, qq/\$(PREFIX).ec.fmd:\$(PREFIX).ec.fq.gz/);
		push(@lines, qq/\t\$(EXE_ROPEBWT2) $opt_rb2 \$< > \$@ 2> \$@.log/, "");

		push(@lines, qq/\$(PREFIX).flt.sub:\$(PREFIX).ec.fmd/);
		push(@lines, qq/\t\$(EXE_FERMI2) occflt -t \$(N_THREADS) -K \$(K_TRIM) \$< > \$@ 2> \$@.log/, "");

		push(@lines, qq/\$(PREFIX).flt.fmd:\$(PREFIX).ec.fmd \$(PREFIX).flt.sub/);
		push(@lines, qq/\t\$(EXE_FERMI2) sub -ct \$(N_THREADS) \$< \$(PREFIX).flt.sub > \$@ 2> \$@.log/, "");
	}

	push(@lines, qq/\$(PREFIX).pre.gz:\$(PREFIX).flt.fmd/);
	push(@lines, qq/\t\$(EXE_FERMI2) assemble -l \$(K_UNITIG) -m \$(K_MERGE) -t \$(N_THREADS) \$< 2> \$@.log | gzip -1 > \$@/, "");

	push(@lines, qq/\$(PREFIX).mag.gz:\$(PREFIX).pre.gz/);
	push(@lines, qq/\t\$(EXE_FERMI2) simplify -CSo \$(K_CLEAN) -m \$(K_MERGE) -T \$(K_UNITIG) \$< 2> \$@.log | gzip -1 > \$@/, "");

	print(join("\n", @lines), "\n");
}

sub utglog {
	die("Usage: fermi2.pl utglog <prefix>\n") if @ARGV == 0;
	while (@ARGV) {
		my $pre = shift(@ARGV);
		my $fh;
		my @a = ($pre, 0, 0, 0, 0);
		open($fh, "$pre.raw.fmd.log") || die;
		while (<$fh>) {
			@a[1,2] = ($1, $2+$3+$4+$5) if /symbol counts.*\((\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+)/;
		}
		close($fh);
		open($fh, "$pre.ec.fq.gz.log") || die;
		while (<$fh>) {
			$a[5] = $1 if /fmc_kmer_stat.*\s(\d+)\s+k/;
		}
		close($fh);
		open($fh, "$pre.ec.fmd.log") || die;
		while (<$fh>) {
			@a[3,4] = ($1, $2+$3+$4+$5) if /symbol counts.*\((\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+)/;
		}
		close($fh);
		if (open($fh, "$pre.mag.gz.log")) {
			while (<$fh>) {
				if (/average read distance ([\d\.]+)/) {
					$a[6] = $1;
				} elsif (/approximate genome size: (\d+)/) {
					$a[7] = $1;
				}
			}
			close($fh);
		} else { $a[6] = 0; $a[7] = 0; }
		print(join("\t", @a), "\n");
	}
}

sub which
{
	my $file = shift;
	my $path = (@_)? shift : $ENV{PATH};
	return if (!defined($path));
	foreach my $x (split(":", $path)) {
		$x =~ s/\/$//;
		return "$x/$file" if (-x "$x/$file");
	}
	return;
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
    } elsif (defined($dirname) && (-x "$dirname/$progname")) {
        return "$dirname/$progname";
    } elsif (-x "./$progname") {
        return "./$progname";
    } elsif (($tmp = &which($progname))) {
        return $tmp;
    } else {
        return;
    }
}

sub dirname {
	my $prog = shift;
	return '.' unless ($prog =~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}
