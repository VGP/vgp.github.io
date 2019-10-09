#!/usr/bin/perl

use strict;
use Time::Local;
#use YAML::Tiny;

#  If set to 1, do not download any genomic_data for coverage estimation.
my $SKIP_RAW = 0;
my $SKIP_ASM = 0;

#  Thresholds for declaring contigs and scaffolds good or bad.
my $goodCTG = 1000000;
my $goodSCF = 10000000;

#  In memory cache of files in genomeark.
#my @genomeArkDates;
#my @genomeArkTimes;
my @genomeArkEpoch;
my @genomeArkSizes;
my @genomeArkFiles;
my $genomeArkLength;

#
#  Discover species.
#

sub discoverDir ($) {
    my $dir = shift @_;
    my @speciesList;

    open(LS, "ls $dir/*.* |");
    while (<LS>) {
        chomp;

        if (m!genomeark/(\w+).md$!)   { push@speciesList, $1; }
        if (m!species/(\w+).yaml$!)   { push@speciesList, $1; }
    }
    close(LS);

    return(@speciesList);
}

sub discover (@) {
    my @species;

    if (scalar(@_) > 0) {
        @species = @_;
    }

    else {
        @species = (discoverDir("vgp-metadata/species"));
        #@species = (discoverDir("../_genomeark"));
    }

    print STDERR "Found ", scalar(@species), " species.\n";

    return(@species);
}


sub loadGenomeArk () {

    #  If no vgp-metadata directory, clone it.

    if (! -e "vgp-metadata") {
        print STDERR "FETCHING METADATA.\n";
        system("git clone git\@github.com:VGP/vgp-metadata.git");
    }

    #  If no genomeark.ls file list, fetch it AND update metadata.

    if (! -e "genomeark.ls") {
        print STDERR "FETCHING AWS FILE LIST.\n";
        system("aws --no-sign-request s3 ls --recursive s3://genomeark/ > genomeark.ls");

        print STDERR "UPDATING METADATA.\n";
        system("cd vgp-metadata ; git fetch ; git merge");
    }

    #  Fail if either doesn't exist.

    die "ERROR: 'genomeark.ls' doesn't exist, can't update.\n"   if (! -e "genomeark.ls");
    die "ERROR: 'vgp-metadata' doesn't exist, can't update.\n"   if (! -e "vgp-metadata");

    #  Pull in all the good bits from the file list.

    open(LS, "< genomeark.ls");
    while (<LS>) {
        chomp;

        my ($filedate, $filetime, $filesize, $filename, $filesecs) = split '\s+', $_;

        my @fileComps   = split '/', $filename;
        my $speciesName = $fileComps[1];
        my $asmName     = $fileComps[3];
        my $seconds     = 0;

        next   if ($filename =~ m!intermediate!);
        next   if ($filename =~ m!transcriptomic_data!);

        if ("$filedate $filetime" =~ m/(\d\d\d\d)-(\d\d)-(\d\d)\s+(\d\d):(\d\d):(\d\d)/) {
            my ($yr, $mo, $da, $hr, $mn, $sc) = ($1, $2, $3, $4, $5, $6);

            $filesecs = timelocal($sc, $mn, $hr, $da, $mo-1, $yr); 
        } else {
            die "Failed to parse date ('$filedate') and time ('$filetime') for file '$filename'\n";
        }

        #push @genomeArkDates, $filedate;
        #push @genomeArkTimes, $filetime;
        push @genomeArkEpoch, $filesecs;
        push @genomeArkSizes, $filesize;
        push @genomeArkFiles, $filename;


        $genomeArkLength++;
    }

    #  Return the time since epoch for the last change to genomeark.ls.

    return((stat("genomeark.ls"))[9]);
}



sub loadAssemblyStatus () {

    my %asmToShow;    #  Map species_name to assembly_name.
    my %asmDate;      #  Map species_name to time_local\0assembly_name.

    #  Scan the genomeark files to figure out what assemblies exist for each
    #  species.  We'll later pick one to use, or use the user-supplied assembly.

    for (my $ii=0; $ii<$genomeArkLength; $ii++) {
        #my $filedate = $genomeArkDates[$ii];
        #my $filetime = $genomeArkTimes[$ii];
        my $filesecs = $genomeArkEpoch[$ii];
        my $filesize = $genomeArkSizes[$ii];
        my $filename = $genomeArkFiles[$ii];

        my @fileComps   = split '/', $filename;
        my $speciesName = $fileComps[1];
        my $asmName     = $fileComps[3];

        next   if ($filename =~ m!genomic_data!);
        next   if ($filename !~ m!fasta.gz!);

        next   if (($filename !~ m/pri.....\d+.fasta/) &&
                   ($filename !~ m/alt.....\d+.fasta/) &&
                   ($filename !~ m/pat.....\d+.fasta/) &&
                   ($filename !~ m/mat.....\d+.fasta/));

        if ($asmName eq "assembly_curated") {   #  If curated exists, always use it,
            $filesecs = 9999999999;             #  even if it is older than some other assembly.
        }

        if ($asmDate{$speciesName} < $filesecs) {
            $asmDate{"$speciesName"} = "$filesecs\0$asmName";
        }
    }

    print STDERR "\n";
    print STDERR "Latest assembly set by file date.\n";
    print STDERR "\n";

    foreach my $speciesName (keys %asmDate) {
        my ($t, $a) = split '\0', $asmDate{$speciesName};

        printf STDERR "%-30s -> %-30s ($t)\n", $speciesName, $a;

        $asmToShow{$speciesName} = $a;
    }

    print STDERR "\n";
    print STDERR "Latest assembly set by user.\n";
    print STDERR "\n";

    open(A, "< assembly_status") or die "Failed to open 'assembly_status' for reading: $!\n";
    while (<A>) {
        chomp;

        next   if (m/^#/);

        if (m/^(\w+_\w+)\s+(a.*)$/) {
            if ($asmToShow{$1} ne $2) {
                printf STDERR "%-30s -> %-30s (previously %s)\n", $1, $2, $asmToShow{$1};
            } else {
                printf STDERR "%-30s -> %-30s\n", $1, $2;
            }

            $asmToShow{$1} = $2;
        }

        else {
            die "Failed to parse 'assembly_status' line '$_'\n";
        }
    }
    close(A);

    return(%asmToShow);
}


#
#  Load the existing .md page for a species.
#

sub loadMeta ($$) {
    my $species = shift @_;
    my $meta    = shift @_;

    my  @keys;
    my  @lvls;

    my  $lvl = 0;

    undef %$meta;

    open(MD, "< vgp-metadata/species/$species.yaml");
    while (<MD>) {
        chomp;

        if (m/^(\s*)(\S*):\s*(.*)$/) {
            my $indent = $1;
            my $key    = $2;   $key   =~ s/^\s+//;  $key   =~ s/\s+$//;
            my $value  = $3;   $value =~ s/^\s+//;  $value =~ s/\s+$//;

            my $len    = length($indent);

            #print STDERR "\n";
            #print STDERR "WORK $len $lvl $key -> $value\n";

            if      ($len  < $lvl) {
                while ($len < $lvl) {
                    #print STDERR "pop     ", $keys[-1], " len=$len lvl=$lvl\n";
                    $lvl -= $lvls[-1];
                    pop @keys;
                    pop @lvls;
                }
            }

            if ($len == $lvl) {
                #print STDERR "replace ", $keys[-1], "\n";
                pop @keys;
                push @keys, $key;

            } elsif ($len >  $lvl) {
                #print STDERR "append  ", $keys[-1], "\n";
                push @keys, $key;
                push @lvls, $len - $lvl;

                $lvl = $len;
            }

            $key = join '.', @keys;

            #print STDERR "$key: $value\n";

            $$meta{$key} = $value;
        }
    }

    die "No meta{species.name} found?\n"  if ($$meta{"species.name"} eq "");

    my @n = split '\s+', $$meta{"species.name"};

    return("$n[0]_$n[1]");      #  Return the proper Species_name for this species.
}



sub loadData ($$) {
    my $species = shift @_;
    my $data    = shift @_;

    my  @keys;
    my  @lvls;

    my  $lvl = 0;

    undef %$data;

    my $key;

    open(MD, "< ../_genomeark/$species.md");
    while (<MD>) {
        chomp;

        #  Match the YAML delimiters.
        if    (m!^---$!) {
            $key    = undef;
        }

        #  Match any multi-line pre-formatted parameters.
        elsif (m!(\w+):\s*\|$!) {
            $key = $1;
            $$data{$key} = "|\n";
        }

        #  Add data to multi-line parameters.
        elsif (defined($key) && (m!^\s\s\S!)) {
            $$data{$key} .= "$_\n";
        }

        #  Match any key-value pairs.
        elsif (m!(\w+):\s*(.*)$!) {
            $key      = undef;
            $$data{$1} = $2;
        }

        #  Match any extra stuff in the file.
        else {
            $$data{":"} .= "$_\n";
        }
    }
    close(MD);
}



sub saveData ($$) {
    my $species = shift @_;
    my $data    = shift @_;

    my @keys = keys %$data;

    @keys = sort @keys;

    open(MD, "> ../_genomeark/$species.md");

    print MD "---\n";
    foreach my $key (@keys) {
        next if ($key eq ":");

        chomp $$data{$key};   #  Remove the newline on multi-line data.

        if ($$data{$key} ne "") {             #  Just because the original didn't
            print MD "$key: $$data{$key}\n";  #  include a space when the data
        } else {                              #  was null.
            #print MD "$key:\n";
        }
    }
    print MD "---\n";
    print MD $$data{":"};

    close(MD);
}



sub printData ($$$) {
    my $species = shift @_;
    my $data    = shift @_;

    my @keys = keys %$data;

    @keys = sort @keys;

    foreach my $key (@keys) {
        next if ($key eq ":");

        chomp $$data{$key};   #  Remove the newline on multi-line data.

        if ($$data{$key} ne "") {
            print "$key: $$data{$key}\n";
        } else {
            print "$key:\n";
        }
    }
    print $$data{":"};
}



sub prettifyBases ($) {
    my $size = shift @_;

    if ($size eq "-")   { return("$size"); }   #  Empty row in N50 table.
    if ($size == 0)     { return(undef);   }   #  Unset genome size.

    #  Shows up to "### bp"
    if      ($size < 1000) {
        return("$size  bp");
    }

    #  Shows up to "##.## Kbp"
    elsif ($size < 100000) {
        return(sprintf("%.2f Kbp", $size / 1000));
    }

    #  Shows up to "500.## Mbp"
    elsif ($size < 500000000) {
        return(sprintf("%.2f Mbp", $size / 1000000));
    }

    #  Everything else.
    else {
        return(sprintf("%.2f Gbp", $size / 1000000000));
    }
}



sub generateAssemblySummary ($$$) {
    my $filename   = shift @_;
    my $type       = shift @_;
    my $genomeSize = shift @_;
    my $split      = ($type eq "ctg") ? "-split-n" : "";


    if (-z "$filename.$type.summary") {
        unlink "$filename.$type.summary";
    }

    if (! -e "downloads") {
        system("mkdir -p downloads");
    }

    if ((! -e "downloads/$filename.gz") &&
        (! -e "$filename.$type.summary")) {
        print STDERR "Fetch s3://genomeark/$filename.gz\n";

        system("aws --no-sign-request s3 cp s3://genomeark/$filename.gz downloads/$filename.gz");
    }

    if ((! -e "downloads/$filename.gz") &&
        (! -e "$filename.$type.summary")) {
        print STDERR "FAILED TO FETCH '$filename'\n";
        exit(1);
    }

    if (! -e "$filename.$type.summary") {
        my $cmd;

        $cmd  = "sequence summarize $split -1x";
        $cmd .= " -size $genomeSize"                                  if ($genomeSize > 0);
        $cmd .= " downloads/$filename.gz > $filename.$type.summary";

        print STDERR "SUMMARIZING '$filename' with genome_size $genomeSize and options '$split'\n";
        print STDERR "  $cmd\n";
        system("mkdir -p $filename");
        system("rmdir    $filename");
        system($cmd);
    }

    if (! -e "$filename.$type.summary") {
        print STDERR "FAILED TO FIND SIZES '$filename'\n";
        exit(1);
    }

    system("rm -f downloads/$filename.gz");
}



sub loadAssemblySummary ($$$$$$) {
    my $filename   = shift @_;
    my $ctgNG      = shift @_;
    my $ctgLG      = shift @_;
    my $ctgLEN     = shift @_;
    my $ctgCOV     = shift @_;
    my $genomeSize = shift @_;
    my $n50        = undef;

    open(SU, "< $filename");
    while (<SU>) {
        s/^\s+//;
        s/\s+$//;
        next if ($_ eq "");
        my ($a, $b) = split '||', $_;
        $a =~ s/^\s+//;
        $a =~ s/\s+$//;

        #  Usually:
        #    a = NG
        #    b = length
        #    c = index
        #    d = sum
        #  except on the last line where the is no length field.
        #
        my ($a, $b, $c, $d) = split '\s+', $_;

        #  If zero, the column was empty, or not a number.
        if ($a == 0) {
            next;
        }

        #  Save the n50 value before we prettify it.
        if (($a == 50) && ($b ne "-")) {
            $n50 = $b;
        }

        #  Add the values to the table, treating the summary line as a footer.
        if ($a =~ m!x!) {
            $a =~ s/^0+//;      #  Strip off leading zeros from "000.851x"
            $a =~ s/^\./0./;    #  But then add back one to get "0.851x"
            $c = prettifyBases($c);

            push @$ctgNG,  $a;
            push @$ctgLG,  $b;
            push @$ctgLEN, $c;
            push @$ctgCOV, $c;

            #$n50table .= "  </tbody>\n";
            #$n50table .= "  <tfoot>\n";
            #$n50table .= "  <tr><th>$a</th><th>$b</th><th></th><th>$c</th></tr>\n";
            #$n50table .= "  </tfoot>\n";
            last;
        }

        else {
            $a =~ s/^0+//;
            #$b = ($b eq "-") ? "-" : sprintf("%.1f", $b / 1000000);
            $b = prettifyBases($b);
            $d = (($genomeSize == 0) || ($b eq "-")) ? "-" : sprintf("%.2f", $d / $genomeSize);

            push @$ctgNG,  $a;
            push @$ctgLG,  $c;
            push @$ctgLEN, $b;
            push @$ctgCOV, $d;

            #my $scfcolor = "";
            #my $ctgcolor = "";

            #$n50table .="  <tr style=\"background-color:#eecccc;\"><td>$a</td><td>$c</td><td><b>$b</b></td><td>$d</td>\n";
            #$n50table .="  <tr style=\"background-color:#cceecc;\"><td>$a</td><td>$c</td><td><b>$b</b></td><td>$d</td>\n";

            #if ($a == 50) {
            #    #($n50 < 5000000)) {
            #    #($n50 >= 5000000)) {
            #}

            #$n50table .= "  <tr><td>$a</td>";
            #$n50table .= "<td>$c</td><td>$b</td><td>$d</td>";
            #$n50table .= "<td>$c</td><td>$b</td><td>$d</td>";
            #$n50table .= "</tr>";

        }
    }
    close(SU);

    return($n50);
}



sub generateAssemblySummaryHTML ($$$) {
    my $prialt     = shift @_;
    my $filename   = shift @_;
    my $genomeSize = shift @_;
    my $n50table   = "";

    $filename =~ s/.gz//;

    generateAssemblySummary($filename, "ctg", $genomeSize);
    generateAssemblySummary($filename, "scf", $genomeSize);

    $n50table .= "|\n";
    $n50table .= "  <table class=\"sequence-sizes-table\">\n";
    #$n50table .= "  <colgroup>\n";
    #$n50table .= "  <col width=\"10%\" />\n";
    #$n50table .= "  <col width=\"10%\" />\n";
    #$n50table .= "  <col width=\"40%\" />\n";
    #$n50table .= "  <col width=\"40%\" />\n";
    #$n50table .= "  </colgroup>\n";
    $n50table .= "  <thead>\n";

    $n50table .= "  <tr>\n";
    $n50table .= "  <th></th>\n";
    $n50table .= "  <th colspan=2 align=center>Contigs</th>\n";
    $n50table .= "  <th colspan=2 align=center>Scaffolds</th>\n";
    $n50table .= "  </tr>\n";

    $n50table .= "  <tr>\n";
    $n50table .= "  <th>NG</th>\n";
    $n50table .= "  <th>LG</th>\n";
    $n50table .= "  <th>Len</th>\n";
    $n50table .= "  <th>LG</th>\n";
    $n50table .= "  <th>Len</th>\n";
    $n50table .= "  </tr>\n";
    $n50table .= "  </thead>\n";
    $n50table .= "  <tbody>\n";

    #  Strip out the N50 chart

    my (@ctgNG, @ctgLG, @ctgLEN, @ctgCOV);
    my (@scfNG, @scfLG, @scfLEN, @scfCOV);

    my $ctgn50 = loadAssemblySummary("$filename.ctg.summary", \@ctgNG, \@ctgLG, \@ctgLEN, \@ctgCOV, $genomeSize);
    my $scfn50 = loadAssemblySummary("$filename.scf.summary", \@scfNG, \@scfLG, \@scfLEN, \@scfCOV, $genomeSize);

    #  Ten rows of actual data.

    for (my $ii=0; $ii<10; $ii++) {
        my $ctgcolor = "";
        my $scfcolor = "";

        #  For non-alternate assemblies (primary, maternal, paternal) highlight the n-50 line.

        if (($ii == 4) && ($prialt ne "alt")) {
            $ctgcolor = ($ctgn50 <  $goodCTG) ? " style=\"background-color:#ff8888;\"" : " style=\"background-color:#88ff88;\"";
            $scfcolor = ($scfn50 <  $goodSCF) ? " style=\"background-color:#ff8888;\"" : " style=\"background-color:#88ff88;\"";
        }

        if ($ii == 4) {
            $n50table .= "  <tr style=\"background-color:#cccccc;\">";
        } else {
            $n50table .= "  <tr>";
        }

        $n50table .= "<td> $ctgNG[$ii] </td>";
        $n50table .= "<td> $ctgLG[$ii] </td><td$ctgcolor> $ctgLEN[$ii] </td>";  #<td> $ctgCOV[$ii] </td>
        $n50table .= "<td> $scfLG[$ii] </td><td$scfcolor> $scfLEN[$ii] </td>";  #<td> $scfCOV[$ii] </td>
        $n50table .= "</tr>";
    }

    #  And one row of summary.

    $n50table .= "  </tbody>\n";
    $n50table .= "  <tfoot>\n";
    $n50table .= "  <tr>";
    $n50table .= "<th> $ctgNG[10] </th>";
    $n50table .= "<th> $ctgLG[10] </th><th> $ctgLEN[10] </th>";
    $n50table .= "<th> $scfLG[10] </th><th> $scfLEN[10] </th>";
    $n50table .= "</tr>\n";
    $n50table .= "  </tfoot>\n";

    $n50table .= "  </table>\n";

    return($ctgn50, $scfn50, $n50table);
}




sub processData ($$$$$) {
    my $filesize = shift @_;
    my $filename = shift @_;
    my $seqFiles = shift @_;
    my $seqBytes = shift @_;
    my $seqIndiv = shift @_;

    my $sTag = "";
    my $iTag = "";

    if ($filename =~ m!species/(.*)/(.*)/genomic_data!) {
        $sTag = $1;
        $iTag = $2;
    } else {
        die "failed to parse species name and individual from '$filename'\n";
    }

    #  Based on directory and filename, count stuff.

    if    ($filename =~ m!/genomic_data/10x/!) {
        $$seqIndiv{"10x"} .= "$sTag/$iTag\0";

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"10x"} += 1;
            $$seqBytes{"10x"} += $filesize;
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/arima/!) {
        $$seqIndiv{"arima"} .= "$sTag/$iTag\0";

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"arima"} += 1;
            $$seqBytes{"arima"} += $filesize;
        } elsif ($filename =~ m/re_bases.txt/) {
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/bionano/!) {
        $$seqIndiv{"bionano"} .= "$sTag/$iTag\0";

        if      ($filename =~ m/cmap/) {
            $$seqFiles{"bionano"} += 1;
            $$seqBytes{"bionano"} += $filesize;
        } elsif ($filename =~ m/bnx.gz/) {
            $$seqBytes{"bionano"} += $filesize;
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/dovetail/!) {
        $$seqIndiv{"dovetail"} .= "$sTag/$iTag\0";

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"dovetail"} += 1;
            $$seqBytes{"dovetail"} += $filesize;
        } elsif ($filename =~ m/re_bases.txt/) {
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/illumina/!) {
        $$seqIndiv{"illumina"} .= "$sTag/$iTag\0";

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"illumina"} += 1;
            $$seqBytes{"illumina"} += $filesize;
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/pacbio/!) {

        if    ($filename =~ m/scraps.bam/) {
            $$seqFiles{"pbscraps"} += 1;
            $$seqBytes{"pbscraps"} += $filesize;
        }
        elsif ($filename =~ m/scraps.bam.pbi/) {
            #$$seqFiles{"pbscraps"} += 1;
            $$seqBytes{"pbscraps"} += $filesize;
        }
        elsif ($filename =~ m/subreads.bam/) {
            $$seqIndiv{"pbsubreads"} .= "$sTag/$iTag\0";
            $$seqFiles{"pbsubreads"} += 1;
            $$seqBytes{"pbsubreads"} += $filesize;
        }
        elsif ($filename =~ m/subreads.bam.pbi/) {
            #$$seqFiles{"pbsubreads"} += 1;
            $$seqBytes{"pbsubreads"} += $filesize;
        }
        elsif ($filename =~ m/ccs.bam/) {
            $$seqIndiv{"pbhifi"} .= "$sTag/$iTag\0";
            $$seqFiles{"pbhifi"} += 1;
            $$seqBytes{"pbhifi"} += $filesize;
        }
        elsif ($filename =~ m/ccs.bam.pbi/) {
            #$$seqFiles{"pbhifi"} += 1;
            $$seqBytes{"pbhifi"} += $filesize;
        }
        elsif ($filename =~ m/ccs.bam.bai/) {
            #$$seqFiles{"pbhifi"} += 1;
            $$seqBytes{"pbhifi"} += $filesize;
        }
        else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/phase/!) {
        $$seqIndiv{"phase"} .= "$sTag/$iTag\0";

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"phase"} += 1;
            $$seqBytes{"phase"} += $filesize;
        } elsif ($filename =~ m/re_bases.txt/) {
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/!) {
        print STDERR "UNKNOWN genomic_data $_\n";
    }
}



sub processAssembly ($$$) {
    my $filesize = shift @_;
    my $filename = shift @_;
    my $data     = shift @_;

    print "PROCESS $filename of size $filesize\n";

    my ($sName, $aLabel, $sTag, $sNum, $prialt, $date) = undef;

    if    ($filename =~ m!species/(.*)/.*/(assembly_.+)/(.......)(\d).(\w+).\w+.(\d\d\d\d)(\d\d)(\d\d).fasta.gz!i) {
        $sName   = $1;
        $aLabel  = $2;
        $sTag    = $3;
        $sNum    = $4;
        $prialt  = $5;
        $date    = "$6-$7-$8";
    }

    #                 species/Gopherus_evgoodei/rGopEvg1/assembly_mt_milan/rGopEvg1.MT.20190310.fasta.gz
    elsif ($filename =~ m!species/(.*)/.*/(assembly_.+)/(.......)(\d).MT.(\d\d\d\d)(\d\d)(\d\d).fasta.gz!i) {
        print STDERR "MITO $1 $2 $3 $4 $6-$7-$8\n";
        $sName   = $1;
        $aLabel  = $2;
        $sTag    = $3;
        $sNum    = $4;
        $prialt  = "mito";
        $date    = "$5-$6-$7";
    }

    else {
        print STDERR "  Nothig here.\n";
        return;
    }

    #  If more than 4, style.scss and _layouts/genomeark.html need updating.
    die "Too many assemblies; update style.scss and _layouts/genomeark.html.\n"   if ($sNum > 4);

    #  If there isn't a date in our database, set it.

    if (!exists($$data{"${prialt}${sNum}date"})) {
        print "  ADD  ${prialt}${sNum}date = $date\n";
        $$data{"${prialt}${sNum}date"}      = $date;
    }

    #  If the date of this file is older than the date of our data, skip it.

    if ($date lt $$data{"${prialt}${sNum}date"}) {
        print "  SKIP ${prialt}${sNum}date = $date\n";
        next;
    }

    #  If the date of this file is newer than the date of our data, reset the date of our data
    #  and wipe out all the old data.

    if ($$data{"${prialt}${sNum}date"} lt $date) {
        print "  RSET ${prialt}${sNum}date = $date\n";

        $$data{"${prialt}${sNum}date"} = $date;
    }

    #  Update everything for this file.

    #print "  add prialt $prialt sNum ${sNum} date $date for $filename\n";

    #$$data{"${prialt}${sNum}"}         = "${aLabel}${sNum}";

    $$data{"${prialt}${sNum}version"}  = "${aLabel}";

    if      ($filesize < 1024 * 1024) {
        $$data{"${prialt}${sNum}filesize"} = sprintf("%.0f KB", $filesize / 1024);
    } elsif ($filesize < 1024 * 1024 * 1024) {
        $$data{"${prialt}${sNum}filesize"} = sprintf("%.0f MB", $filesize / 1024 / 1024);
    } elsif ($filesize < 1024 * 1024 * 1024 * 1024) {
        $$data{"${prialt}${sNum}filesize"} = sprintf("%.0f GB", $filesize / 1024 / 1024 / 1024);
    } else {
        $$data{"${prialt}${sNum}filesize"} = sprintf("%.0f TB", $filesize / 1024 / 1024 / 1024 / 1024);
    }

    $$data{"${prialt}${sNum}seq"}      = "https://s3.amazonaws.com/genomeark/$filename";

    ($$data{"${prialt}${sNum}n50ctg"},
     $$data{"${prialt}${sNum}n50scf"},
     $$data{"${prialt}${sNum}sizes"})  = generateAssemblySummaryHTML($prialt, $filename, $$data{"genome_size"});

    #  Update the assembly status based on the primary n50 and/or curation status.

    #print STDERR "$prialt -- " . $$data{"${prialt}${sNum}n50ctg"} . " -- " . $$data{"${prialt}${sNum}n50scf"} . "\n";

    if (($prialt eq "pri") ||
        ($prialt eq "mat") ||
        ($prialt eq "pat")) {
        if (($$data{"${prialt}${sNum}n50ctg"} < $goodCTG) ||
            ($$data{"${prialt}${sNum}n50scf"} < $goodSCF)) {
            $$data{"assembly_status"} = "low-quality-draft";
        }

        elsif ($aLabel !~ m/curated/) {
            $$data{"assembly_status"} = "high-quality-draft";
        }

        else {
            $$data{"assembly_status"} = "curated";
        }
    }
}



sub downloadAndSummarize ($$$) {
    my $name  = shift @_;
    my $size  = shift @_;
    my $file  = shift @_;
    my $bases = 0;

    if (-z "$name.summary") {
        unlink "$name.summary";
    }

    if (! -e "downloads") {
        system("mkdir -p downloads");
    }

    if ((! -e "downloads/$name") &&
        (! -e "$name.summary")) {
        printf STDERR "FETCH file #%4d size %6.3f GB '%s'\n", $file, $size / 1024 / 1024 / 1024, $name;
        printf STDERR " -> downloads/$name\n";
        #printf STDERR "  aws --no-sign-request s3 cp s3://genomeark/$name downloads/$name\n";
        system("aws --no-sign-request s3 cp s3://genomeark/$name downloads/$name")  if ($SKIP_RAW == 0);
    }

    #  Make the directory where we'll store summaries.

    #  If a bam, convert to fastq then summarize.

    if ($name =~ m/bam$/) {
        if ((  -e "downloads/$name") &&
            (! -e "$name.summary")) {
            printf STDERR "EXTRACT bases\n";
            system("samtools bam2fq downloads/$name > downloads/$name.fastq");
        }
        if ((  -e "downloads/$name.fastq") &&
            (! -e "$name.summary")) {
            printf STDERR "SUMMARIZE $name.summary\n";
            #printf STDERR "  sequence summarize downloads/$name.fastq > $name.summary\n";
            system("mkdir -p $name");
            system("rmdir    $name");
            system("sequence summarize downloads/$name.fastq > $name.summary");
        }
    }

    #  Otherwise, summarize directly.

    else {
        if ((  -e "downloads/$name") &&
            (! -e "$name.summary")) {
            printf STDERR "SUMMARIZE $name.summary\n";
            #printf STDERR "  sequence summarize downloads/$name > $name.summary\n";
            system("mkdir -p $name");
            system("rmdir    $name");
            system("sequence summarize downloads/$name > $name.summary");
        }
    }

    #  Parse the summary to find the number of bases in the dataset.

    if (  -e "$name.summary") {
        open(ST, "< $name.summary");
        while (<ST>) {
            if (m/^G=(\d+)\s/) {
                $bases = $1;
                last;
            }
        }
        close(ST);
    }

    if (($bases == 0) && ($SKIP_RAW == 1)) {
        return($size);
    }

    return($bases);
}



sub estimateRawDataScaling ($$) {
    my $name    = shift @_;
    my $type    = shift @_;
    my @files;

    for (my $ii=0; $ii<$genomeArkLength; $ii++) {
        #my $filedate = $genomeArkDates[$ii];
        #my $filetime = $genomeArkTimes[$ii];
        my $filesize = $genomeArkSizes[$ii];
        my $filename = $genomeArkFiles[$ii];

        next   if ($filename !~ m!$name!);

        if    (($type eq "10x") && ($filename =~ m!genomic_data/10x/.*q.gz$!)) {
            push @files, "$filesize\0$filename";
        }

        if    (($type eq "arima") && ($filename =~ m!genomic_data/arima/.*q.gz$!)) {
            push @files, "$filesize\0$filename";
        }

        if    (($type eq "dovetail") && ($filename =~ m!genomic_data/dovetail/.*q.gz$!)) {
            push @files, "$filesize\0$filename";
        }

        if    (($type eq "illumina") && ($filename =~ m!genomic_data/illumina/.*q.gz$!)) {
            push @files, "$filesize\0$filename";
        }

        if    (($type eq "pbsubreads") && ($filename =~ m!genomic_data/pacbio/.*subreads.*bam$!)) {
            push @files, "$filesize\0$filename";
        }

        if    (($type eq "pbhifi") && ($filename =~ m!genomic_data/pacbio/.*ccs.*bam$!)) {
            push @files, "$filesize\0$filename";
        }

        if    (($type eq "phase") && ($filename =~ m!genomic_data/phase/.*q.gz$!)) {
            push @files, "$filesize\0$filename";
        }
    }

    @files = sort { $a <=> $b } @files;

    my $nFiles = scalar(@files);

    if ($nFiles == 0) {
        return(1.0);
    }

    my $f1 = int(1 * $nFiles / 4);
    my $f2 = int(2 * $nFiles / 4);
    my $f3 = int(3 * $nFiles / 4);

    my ($size1, $name1, $bases1, $seqs1) = split '\0', @files[$f1];
    my ($size2, $name2, $bases2, $seqs2) = split '\0', @files[$f2];
    my ($size3, $name3, $bases3, $seqs3) = split '\0', @files[$f3];

    $bases1 = downloadAndSummarize($name1, $size1, $f1);
    $bases2 = downloadAndSummarize($name2, $size2, $f2);
    $bases3 = downloadAndSummarize($name3, $size3, $f3);

    if (($bases1 == 0) ||
        ($bases2 == 0) ||
        ($bases3 == 0)) {
        print STDERR "FAILED TO ESTIMATE SIZES.\n";
        print STDERR "  1 - $bases1 - $name1\n";
        print STDERR "  2 - $bases2 - $name2\n";
        print STDERR "  3 - $bases3 - $name3\n";
        die "\n";
    }

    #print STDERR "\n";
    #print STDERR "SCALING 1:  size $size1 -- bases $bases1 -- $name1\n";
    #print STDERR "SCALING 2:  size $size2 -- bases $bases2 -- $name2\n";
    #print STDERR "SCALING 3:  size $size3 -- bases $bases3 -- $name3\n";

    my $scaling = ($bases1 + $bases2 + $bases3) / ($size1 + $size2 + $size3);

    $scaling = int($scaling * 10000) / 10000;  #  Prevent churn by limiting precision of final value.

    return($scaling);
}



sub computeBionanoBases ($) {
    my $name    = shift @_;
    my @files;

    for (my $ii=0; $ii<$genomeArkLength; $ii++) {
        #my $filedate = $genomeArkDates[$ii];
        #my $filetime = $genomeArkTimes[$ii];
        my $filesize = $genomeArkSizes[$ii];
        my $filename = $genomeArkFiles[$ii];

        if    ($filename =~ m!$name.*genomic_data/bionano/.*bnx.gz$!) {
            push @files, "$filesize\0$filename";
        }
    }

    #  Download all the cmap files.

    foreach my $f (@files) {
        my ($size, $name) = split '\0', $f;

        if (! -e "downloads") {
            system("mkdir -p downloads");
        }

        if ((! -e "downloads/$name") &&
            (! -e "$name.summary")) {
            printf STDERR "FETCH size %6.3f GB '%s'\n", $size / 1024 / 1024 / 1024, $name;
            system("aws --no-sign-request s3 cp s3://genomeark/$name downloads/$name")  if ($SKIP_RAW == 0);
        }
    }

    #  Parse each bnx to find the molecule sizes

    foreach my $f (@files) {
        my ($size, $name) = split '\0', $f;

        next  if (! -e "downloads/$name");
        next  if (  -e "$name.summary");

        print STDERR "PARSING $name\n";

        my $nMolecules = 0;
        my $sLength    = 0;

        system("mkdir -p $name");
        system("rmdir    $name");

        open(B, "gzip -dc downloads/$name |");
        while (<B>) {
            if (m/^0\s/) {
                my @v = split '\s+', $_;

                $nMolecules += 1;
                $sLength    += $v[2];
            }
        }
        close(B);

        open(S, "> $name.summary") or die "Failed to open '$name.summary' for writing: $!\n";;
        print S "Molecules: $nMolecules\n";
        print S "Length:    $sLength\n";
        close(S);
    }

    #  Read the summaries to find the bases.

    my $bases = 0;

    foreach my $f (@files) {
        my ($size, $name) = split '\0', $f;

        if (-e "$name.summary") {
            open(S, "< $name.summary") or die "Failed to open '$name.summary' for reading: $!\n";
            while (<S>) {
                $bases += $1   if (m/Length:\s+(.*)$/);
            }
            close(S);
        }
    }

    return($bases);
}



#
#  Main
#

my $lastUpdate = loadGenomeArk();

my @speciesList = discover(@ARGV);
my %asmToShow   = loadAssemblyStatus();


#  Rebuild each species.md file.

foreach my $species (@speciesList) {
    my %seqFiles;
    my %seqBytes;
    my %seqIndiv;
    my %meta;
    my %data;

    #
    #  Load metadata and copy the good bits.
    #

    my $name = loadMeta($species, \%meta);
    my $asm;

    $data{"name"}                = $meta{"species.name"};         #  Name with a space:     'Species_name'
    $data{"name_"}               = $name;                         #  Name with underscore:  'Species_name'
    $data{"short_name"}          = $meta{"species.short_name"};

    $data{"common_name"}         = $meta{"species.common_name"};
    $data{"taxon_id"}            = $meta{"species.taxon_id"};

    #$data{"s3"}                  = "s3://genomeark/species/$name";

    $data{"genome_size"}         = $meta{"species.genome_size"};
    $data{"genome_size_display"} = prettifyBases($meta{"species.genome_size"});  #  Updated later, too.
    $data{"genome_size_method"}  = $meta{"species.genome_size_method"};

    $data{"data_status"}         = "none";  #<em style=\"color:red\">no data</em>";
    $data{"assembly_status"}     = "none";  #<em style=\"color:red\">no assembly</em>";

    #$data{"last_raw_data"}       = Do not set here; should only be present if data exists.
    $data{"last_updated"}        = 0;

    #
    #  Find which assembly to use.
    #

    $data{"assembly"} = $asm    = $asmToShow{$name};

    print STDERR "\n";
    print STDERR "$name -- $asm\n";

    for (my $ii=0; $ii<$genomeArkLength; $ii++) {
        #my $filedate = $genomeArkDates[$ii];
        #my $filetime = $genomeArkTimes[$ii];
        my $filesecs = $genomeArkEpoch[$ii];
        my $filesize = $genomeArkSizes[$ii];
        my $filename = $genomeArkFiles[$ii];

        next   if ($filename !~ m!$name!);

        #print "  '$filename'\n";

        if ($data{"last_updated"} < $filesecs) {
            $data{"last_updated"} = $filesecs;
        }

        #  Process raw data.
        if ($filename =~ m!genomic_data!) {
            if ($data{"last_raw_data"} < $filesecs) {     #  If this isn't set, Raw Data shows
                $data{"last_raw_data"} = $filesecs;       #  "No data.".
            }

            processData($filesize, $filename, \%seqFiles, \%seqBytes, \%seqIndiv);
            next;
        }

        #  If this is the assembly_to_show, process it.
        if (($asm ne "") &&
            ($filename =~ m!$asm.*fasta!i)) {
            print "  '$filename' -- GENOME\n";
            processAssembly($filesize, $filename, \%data);
            next;
        }

        #  But always process mito contigs.
        if ($filename =~ m!assembly_mt.*MT.*fasta!i) {
            print "  '$filename' -- MITO\n";
            processAssembly($filesize, $filename, \%data);
            next;
        }
    }

    #
    #  Track down any link to NCBI.
    #

    open(GB, "< genbank.map");
    while (<GB>) {
        chomp;

        my ($gbacc, $gbnam, $gbtyp);

        if ($_ =~ m/^(GCA_\d+.\d+)\s+(\S+)\s+(.*)$/) {
            $gbacc = $1;
            $gbnam = $2;
            $gbtyp = $3;

            if ($gbtyp =~ m/alt/) {
                $gbtyp = "alt";
            } else {
                $gbtyp = "pri";
            }

            #print "$name $data{'short_name'} -> $gbacc $gbnam $gbtyp\n";

            #  Fuzzy match our short_name against ncbi's gbnam.  Fuzzy because
            #  genbank isn't always correct, for example:
            #    fGadMor -- gadMor3.0
            #    mLynCan -- mLynCan4_v1.alt
            #
            if ($data{'short_name'} =~ m/^.(......)$/) {
                my $tag = $1;

                if ($gbnam =~ m/$tag/i) {
                    if ($gbtyp eq "pri") {
                        #print "MATCH PRI $tag -- $gbnam -- $gbacc\n";
                        $data{'genbank_pri'} = $gbacc;
                    } else {
                        #print "MATCH ALT $tag -- $gbnam -- $gbacc\n";
                        $data{'genbank_alt'} = $gbacc;
                    }
                }
            } else {
                die "Malformed short_name '$data{'short_name'}'.\n";
            }

        } else {
            die "Failed to parse genbank.map line '$_'.\n";
        }
    }
    close(GB);

    #
    #  Finalize the genomic_data by adding to %data.
    #

    #printf "\n";
    #printf "%14d gzip bytes in %3d 10x datasets.\n",               $seqBytes{"10x"},          $seqFiles{"10x"}         if ($seqFiles{"10x"} > 0);
    #printf "%14d gzip bytes in %3d Arima datasets.\n",             $seqBytes{"arima"},        $seqFiles{"arima"}       if ($seqFiles{"arima"} > 0);
    #printf "%14d gzip bytes in %3d BioNano datasets.\n",           $seqBytes{"bionano"},      $seqFiles{"bionano"}     if ($seqFiles{"bionano"} > 0);
    #printf "%14d gzip bytes in %3d Dovetail datasets.\n",          $seqBytes{"dovetail"} / 2, $seqFiles{"dovetail"}    if ($seqFiles{"dovetail"} > 0);
    #printf "%14d gzip bytes in %3d Illumina datasets.\n",          $seqBytes{"illumina"} / 2, $seqFiles{"illumina"}    if ($seqFiles{"illumina"} > 0);
    #printf "%14d bam  bytes in %3d PacBio datasets (scraps).\n",   $seqBytes{"pbscraps"},     $seqFiles{"pbscraps"}    if ($seqFiles{"pbscraps"} > 0);
    #printf "%14d bam  bytes in %3d PacBio datasets (subreads).\n", $seqBytes{"pbsubreads"},   $seqFiles{"pbsubreads"}  if ($seqFiles{"pbsubreads"} > 0);
    #printf "%14d bam  bytes in %3d PacBio datasets (hifi).\n",     $seqBytes{"pbhifi"},       $seqFiles{"pbhifi"}      if ($seqFiles{"pbhifi"} > 0);
    #printf "%14d gzip bytes in %3d Phase datasets.\n",             $seqBytes{"phase"} / 2,    $seqFiles{"phase"}       if ($seqFiles{"phase"} > 0);

    $seqFiles{"dovetail"} /= 2;
    $seqFiles{"illumina"} /= 2;
    $seqFiles{"phase"}    /= 2;

    foreach my $type (qw(10x arima bionano dovetail illumina pbscraps pbsubreads pbhifi phase)) {
        if ($seqBytes{$type} > 0) {
            #$data{"data_bytes"}           += $seqBytes{$type};
            $data{"data_${type}_bytes"}    = sprintf("%.3f GB", $seqBytes{$type} / 1024 / 1024 / 1024);
            $data{"data_${type}_coverage"} = "N/A";
            $data{"data_${type}_bases"}    = "unknown";
            $data{"data_${type}_files"}    = $seqFiles{$type};
        }
    }

    if ($data{"data_10x_scale"} == 0) {
        $data{"data_10x_scale"} = estimateRawDataScaling($name, "10x");
    }

    if ($data{"data_arima_scale"} == 0) {
        $data{"data_arima_scale"} = estimateRawDataScaling($name, "arima");
    }

    #  BIONANO is totally different.

    if ($data{"data_dovetail_scale"} == 0) {
        $data{"data_dovetail_scale"} = estimateRawDataScaling($name, "dovetail");
    }

    if ($data{"data_illumina_scale"} == 0) {
        $data{"data_illumina_scale"} = estimateRawDataScaling($name, "illumina");
    }

    if ($data{"data_pbsubreads_scale"} == 0) {
        $data{"data_pbsubreads_scale"} = estimateRawDataScaling($name, "pbsubreads");
    }

    if ($data{"data_pbhifi_scale"} == 0) {
        $data{"data_pbhifi_scale"} = estimateRawDataScaling($name, "pbhifi");
    }

    if ($data{"data_phase_scale"} == 0) {
        $data{"data_phase_scale"} = estimateRawDataScaling($name, "phase");
    }

    #  If no genome size set, default to assembly size.

    if ($data{"genome_size"} == 0)   { $data{"genome_size"} = $data{"pri1length"}; }
    if ($data{"genome_size"} == 0)   { $data{"genome_size"} = $data{"pri2length"}; }
    if ($data{"genome_size"} == 0)   { $data{"genome_size"} = $data{"pri3length"}; }
    if ($data{"genome_size"} == 0)   { $data{"genome_size"} = $data{"pri4length"}; }

    $data{"genome_size_display"} = prettifyBases($data{"genome_size"});

    #  Figure out how much and what types of data exist.

    my $dataPac = 0;
    my $data10x = 0;
    my $dataHIC = 0;
    my $dataBio = 0;

    sub uniquifyStringArray ($) {     #  Ain't perl great!
        my @a = split '\0', $_[0];
        my %a;

        foreach my $a (@a) {
            $a{$a}++;
        }

        return(sort keys %a);
    }


    if (($seqBytes{"10x"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"10x"})) {
            $data{"data_10x_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/10x/ .<br>";
        }

        if ($data{"genome_size"} > 0) {
            $data{"data_10x_coverage"}        = sprintf("%.2fx", $seqBytes{"10x"} * $data{"data_10x_scale"} / $data{"genome_size"});
            $data{"data_10x_bases"}           = prettifyBases($seqBytes{"10x"} * $data{"data_10x_scale"});
        }
        $data10x++;
    }

    if (($seqBytes{"arima"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"arima"})) {
            $data{"data_arima_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/arima/ .<br>";
        }

        if ($data{"genome_size"} > 0) {
            $data{"data_arima_coverage"}      = sprintf("%.2fx", $seqBytes{"arima"} * $data{"data_arima_scale"} / $data{"genome_size"});
            $data{"data_arima_bases"}         = prettifyBases($seqBytes{"arima"} * $data{"data_arima_scale"});
        }
        $dataHIC++;
    }

    if (($seqBytes{"bionano"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"bionano"})) {
            $data{"data_bionano_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/bionano/ .<br>";
        }

        if ($data{"genome_size"} > 0) {
            my $b = computeBionanoBases($name);
            $data{"data_bionano_coverage"}    = sprintf("%.2fx", $b / $data{"genome_size"});
            $data{"data_bionano_bases"}       = prettifyBases($b);
        }
        $dataBio++;
    }

    if (($seqBytes{"dovetail"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"dovetail"})) {
            $data{"data_dovetail_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/dovetail/ .<br>";
        }

        if ($data{"genome_size"} > 0) {
            $data{"data_dovetail_coverage"}   = sprintf("%.2fx", $seqBytes{"dovetail"} * $data{"data_dovetail_scale"} / $data{"genome_size"});
            $data{"data_dovetail_bases"}      = prettifyBases($seqBytes{"dovetail"} * $data{"data_dovetail_scale"});
        }
        $dataHIC++;
    }

    if (($seqBytes{"illumina"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"illumina"})) {
            $data{"data_illumina_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/illumina/ .<br>";
        }

        if ($data{"genome_size"} > 0) {
            $data{"data_illumina_coverage"}   = sprintf("%.2fx", $seqBytes{"illumina"} * $data{"data_illumina_scale"} / $data{"genome_size"});
            $data{"data_illumina_bases"}      = prettifyBases($seqBytes{"illumina"} * $data{"data_illumina_scale"});
        }
    }

    if (($seqBytes{"pbsubreads"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"pbsubreads"})) {
            $data{"data_pbsubreads_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/pacbio/ . --exclude \"*ccs.bam*\"<br>";
        }

        if ($data{"genome_size"} > 0) {
            $data{"data_pbsubreads_coverage"} = sprintf("%.2fx", $seqBytes{"pbsubreads"} * $data{"data_pbsubreads_scale"} / $data{"genome_size"});
            $data{"data_pbsubreads_bases"}    = prettifyBases($seqBytes{"pbsubreads"} * $data{"data_pbsubreads_scale"});
        }
        $dataPac++;
    }

    if (($seqBytes{"pbhifi"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"pbhifi"})) {
            $data{"data_pbhifi_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/pacbio/ . --exclude \"*subreads.bam*\"<br>";
        }

        if ($data{"genome_size"} > 0) {
            $data{"data_pbhifi_coverage"} = sprintf("%.2fx", $seqBytes{"pbhifi"} * $data{"data_pbhifi_scale"} / $data{"genome_size"});
            $data{"data_pbhifi_bases"}    = prettifyBases($seqBytes{"pbhifi"} * $data{"data_pbhifi_scale"});
        }
        $dataPac++;
    }

    if (($seqBytes{"phase"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"phase"})) {
            $data{"data_phase_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/phase/ .<br>";
        }

        if ($data{"genome_size"} > 0) {
            $data{"data_phase_coverage"}      = sprintf("%.2fx", $seqBytes{"phase"} * $data{"data_phase_scale"} / $data{"genome_size"});
            $data{"data_phase_bases"}         = prettifyBases($seqBytes{"phase"} * $data{"data_phase_scale"});
        }
        $dataHIC++;
    }

    #  Set the classification.
    #  If no assembly, put it under either "some data" or "all data".
    #  Otherwise, put it under the assembly classification.

    if (($dataPac > 0) ||
        ($data10x > 0) ||
        ($dataHIC > 0) ||
        ($dataBio > 0)) {
        $data{"data_status"} = "some";
    }

    if (($dataPac > 0) &&
        ($data10x > 0) &&
        ($dataHIC > 0) &&
        ($dataBio > 0)) {
        $data{"data_status"} = "all";
    }

    #  Create symlinks to the categories.
    #
    #  Name changes should change:
    #  1)  the label used in this script ("low-quality-draft")
    #  2)  the directory name in the symlink
    #  3)  the actual directories (both "_genomeark-low-quality-draft" and "genomeark-low-quality-draft")
    #  4)  the index.html in genomeark-low-quality-draft/  (both the title and the collection name)
    #  5)  _config.yml (in two places) and index.html in the root directory.

    unlink("../_genomeark-curated-assembly/$name.md");
    unlink("../_genomeark-high-quality-draft-assembly/$name.md");
    unlink("../_genomeark-low-quality-draft-assembly/$name.md");
    unlink("../_genomeark-complete-data/$name.md");
    unlink("../_genomeark-incomplete-data/$name.md");

    if    ($data{"assembly_status"} eq "curated") {
        system("mkdir -p ../_genomeark-curated-assembly");
        symlink("../_genomeark/$name.md", "../_genomeark-curated-assembly/$name.md") or die "Failed to make symlink for curated assembly: $!.\n";
    }

    elsif ($data{"assembly_status"} eq "high-quality-draft") {
        system("mkdir -p ../_genomeark-high-quality-draft-assembly");
        symlink("../_genomeark/$name.md", "../_genomeark-high-quality-draft-assembly/$name.md") or die "Failed to make symlink for high-quality-draft assembly: $!.\n";
    }

    elsif ($data{"assembly_status"} eq "low-quality-draft-assembly") {
        system("mkdir -p ../_genomeark-low-quality-draft");
        symlink("../_genomeark/$name.md", "../_genomeark-low-quality-draft-assembly/$name.md") or die "Failed to make symlink for low-quaity-draft assembly: $!.\n";
    }

    elsif ($data{"data_status"} eq "all") {
        system("mkdir -p ../_genomeark-complete-data");
        symlink("../_genomeark/$name.md", "../_genomeark-complete-data/$name.md") or die "Failed to make symlink for complete data: $!.\n";
    }

    elsif ($data{"data_status"} eq "some") {
        system("mkdir -p ../_genomeark-incomplete-data");
        symlink("../_genomeark/$name.md", "../_genomeark-incomplete-data/$name.md") or die "Failed to make symlink for incomplete data: $!.\n";
    }

    else {
        system("mkdir -p ../_genomeark-incomplete-data");
        symlink("../_genomeark/$name.md", "../_genomeark-incomplete-data/$name.md") or die "Failed to make symlink for incomplete data (catch all): $!.\n";
    }

    #  And reset the classifications to strings we can use in the display.

    $data{"data_status"}     = "<em style=\"color:red\">no data</em>"                          if ($data{"data_status"}) eq "none";
    $data{"data_status"}     = "<em style=\"color:orange\">some data</em>"                     if ($data{"data_status"}) eq "some";
    $data{"data_status"}     = "<em style=\"color:green\">all data</em>"                       if ($data{"data_status"}) eq "all";

    $data{"assembly_status"} = "<em style=\"color:red\">no assembly</em>"                      if ($data{"assembly_status"} eq "none");
    $data{"assembly_status"} = "<em style=\"color:red\">low-quality draft assembly</em>"       if ($data{"assembly_status"} eq "low-quality-draft");
    $data{"assembly_status"} = "<em style=\"color:orange\">high-quality draft assembly</em>"   if ($data{"assembly_status"} eq "high-quality-draft");
    $data{"assembly_status"} = "<em style=\"color:green\">curated assembly</em>"               if ($data{"assembly_status"} eq "curated");

    #  If no date set -- no raw data, no assemblies, no anything -- set it to the
    #  date of this update (genomeark.ls's date).

    if ($data{"last_updated"} == 0) {
        $data{"last_updated"} = $lastUpdate;
    }

    #  Done.  Write the output.

    #printData($species, \%data);
    saveData($species, \%data);
}
