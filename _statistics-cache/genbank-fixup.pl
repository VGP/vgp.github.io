#!/usr/bin/env perl

#  Use genbank-fetch.sh to pull the genbank XML and parse out a few fields.
#  Use genbank-fixup.pl to clean up the result and write genbank.map.

use strict;

open(F, "< genbank.map.raw") or die;
open(O, "> genbank.map") or die;
while (<F>) {
    chomp;

    if ($_ =~ m/^(GCA_\d+.\d+)\s+(\S+)\s+(.*)$/) {
        my $acc = $1;  #  Correct as is.
        my $lab = $2;  #  Minor mods needed.
        my $dsc = $3;  #  A mess.

        my $nam = "XXX";  #  Properly formatted VGP individual name
        my $typ = "XXX";  #  Type of this assembly (pri alt mat pat)

        #  Fix up obvious mistakes in the name.

        if ($lab =~ m/g(adMor.*)/)       { $lab = "fG$1"; }
        if ($lab =~ m/mCalJa(1.2.*)/)    { $lab = "mCalJac$1"; }

        if ($lab =~ m/mRatBN7/)          { next; }
        if ($lab =~ m/UOA_/)             { next; }
        if ($lab =~ m/ZJU/)              { next; }

        if ($lab =~ m/^([a-z][A-Z][a-z][a-z][a-zA-z][A-Za-z][A-Za-z][0-9])/) {
            $nam = $1;  #  matches fAaaAaa or fAaaaAa
        } else {
            die "Failed to match name in '$lab'.\n";
        }

        #  See if the description tells us primary or alternate.
        #  There seem to only be three forms here:
        #    % cut -f 3 genbank.map.full | sort | uniq -c
        #    70 alternate-pseudohaplotype
        #    27 haploid
        #    65 haploid (principal pseudohaplotype of diploid)

        if (($dsc =~ m/princ/) ||
            ($lab =~ m/\.pri/)) {
            $typ = "pri";
        }
        if (($dsc =~ m/alter/) ||
            ($lab =~ m/\.alt/)) {
            $typ = "alt";
        }

        #  See if the name is indicating paternal or maternal

        if (($lab =~ m/\.mat$/) ||
            ($dsc =~ m/maternal/)) {
            $typ = "mat";
        }
        if (($lab =~ m/\.pat$/) ||
            ($dsc =~ m/paternal/)) {
            $typ = "pat";
        }

        #  Hardcode some ugly ones.

        if ($lab eq "mArvNil1.pat.X")     { $typ = "XXX"; }

        if ($lab eq "bTaeGut2.pat.W.v2")  { $typ = "XXX"; }

        if ($lab eq "bTaeGut2.mat.v3")    { $typ = "mat"; }
        if ($lab eq "bTaeGut2pat")        { $typ = "pat"; }

        if ($lab eq "bTaeGut2.pri.v2")    { $typ = "pri"; }
        if ($lab eq "bTaeGut2.p.v1.alt")  { $typ = "alt"; }

        if ($lab eq "eAstRub1.3")   { $typ = "pri"; }
        if ($lab eq "fAstCal1.2")   { $typ = "pri"; }
        if ($lab eq "fBetSpl5.2")   { $typ = "pri"; }
        if ($lab eq "fSalTru1.1")   { $typ = "pri"; }

        #  Emit output or fail.
        if ($typ eq "XXX") {
            print "$acc  $nam  $typ -- $lab -- $dsc -- ignored\n";
        } else {
            print O "$acc\t$nam\t$typ\n";
        }
    }

    else {
        die "Failed to match '$_'\n";
    }
}
close(O);
close(F);
