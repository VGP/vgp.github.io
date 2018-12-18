#!/usr/bin/perl

use strict;
use Time::Local;
#use YAML::Tiny;

#  If set to 1, do not download any genomic_data for coverage estimation.
my $SKIP_RAW = 1;
my $SKIP_ASM = 0;

#  Thresholds for declaring contigs and scaffolds good or bad.
my $goodCTG = 1000000;
my $goodSCF = 10000000;



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


sub loadAssemblyStatus () {
    my %asmToShow;

    open(A, "< assembly_status") or die "Failed to open 'assembly_status' for reading: $!\n";
    while (<A>) {
        chomp;

        if      (m/^(\w+_\w+)\s*$/) {
            $asmToShow{$1} = "";
        }

        elsif (m/^(\w+_\w+)\s+(a.*)$/) {
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
    my $data    = shift @_;

    my  @keys;
    my  @lvls;

    my  $lvl = 0;

    undef %$data;

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

            $$data{$key} = $value;
        }
    }
}


#  Returns the name of the directory this species is in, based on the
#  full species name in the metadata.
sub makeName ($) {
    my $name = shift @_;
    my @n = split '\s+', $name;

    return("$n[0]_$n[1]");

#    if ($name =~ m/(\S+)\s+(\S+)\s+\S+/) {
#        $name = "$1_$2";
#    }
#
#    $name =~ s/\s+/_/g;
#
#    return($name);
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
        print STDERR "SUMMARIZING '$filename' with genome_size $genomeSize and options '$split'\n";
        system("mkdir -p $filename");
        system("rmdir    $filename");
        system("sequence summarize $split -1x -size $genomeSize downloads/$filename.gz > $filename.$type.summary");
    }

    if (! -e "$filename.$type.summary") {
        print STDERR "FAILED TO FIND SIZES '$filename'\n";
        exit(1);
    }
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
            $d = ($b eq "-") ? "-" : sprintf("%.2f", $d / $genomeSize);

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

        if (($ii == 4) && ($prialt eq "pri")) {
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




sub processData ($$$$) {
    my $filesize = shift @_;
    my $filename = shift @_;
    my $seqFiles = shift @_;
    my $seqBytes = shift @_;

    my $sTag = "";   #  Not used here.

    #  Based on directory and filename, count stuff.

    if    ($filename =~ m!/genomic_data/10x/!) {
        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"10x$sTag"} += 1;
            $$seqBytes{"10x$sTag"} += $filesize;
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/arima/!) {
        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"arima$sTag"} += 1;
            $$seqBytes{"arima$sTag"} += $filesize;
        } elsif ($filename =~ m/re_bases.txt/) {
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/bionano/!) {
        if      ($filename =~ m/cmap/) {
            $$seqFiles{"bionano$sTag"} += 1;
            $$seqBytes{"bionano$sTag"} += $filesize;
        } elsif ($filename =~ m/bnx.gz/) {
            $$seqBytes{"bionano$sTag"} += $filesize;
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/dovetail/!) {
        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"dovetail$sTag"} += 1;
            $$seqBytes{"dovetail$sTag"} += $filesize;
        } elsif ($filename =~ m/re_bases.txt/) {
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/illumina/!) {
        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"illumina$sTag"} += 1;
            $$seqBytes{"illumina$sTag"} += $filesize;
        } else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/pacbio/!) {
        if    ($filename =~ m/scraps.bam/) {
            $$seqFiles{"pbscraps$sTag"} += 1;
            $$seqBytes{"pbscraps$sTag"} += $filesize;
        }
        elsif ($filename =~ m/scraps.bam.pbi/) {
            #$$seqFiles{"pbscraps$sTag"} += 1;
            $$seqBytes{"pbscraps$sTag"} += $filesize;
        }
        elsif ($filename =~ m/subreads.bam/) {
            $$seqFiles{"pbsubreads$sTag"} += 1;
            $$seqBytes{"pbsubreads$sTag"} += $filesize;
        }
        elsif ($filename =~ m/subreads.bam.pbi/) {
            #$$seqFiles{"pbsubreads$sTag"} += 1;
            $$seqBytes{"pbsubreads$sTag"} += $filesize;
        }
        else {
            print STDERR "unknown file type in '$_'\n";
        }
    }

    elsif ($filename =~ m!/genomic_data/phase/!) {
        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"phase$sTag"} += 1;
            $$seqBytes{"phase$sTag"} += $filesize;
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

    #print "PROCESS $filename of size $filesize\n";

    my ($aLabel, $sTag, $sNum, $prialt, $date) = undef;

    if    ($filename =~ m!assembly_(.+)/(.......)(\d).(\w+).\w+.(\d\d\d\d)(\d\d)(\d\d).fasta.gz!) {
        $aLabel  = $1;
        $sTag    = $2;
        $sNum    = $3;
        $prialt  = $4;
        $date    = "$5-$6-$7";
    }

    elsif ($filename =~ m!assembly_(,+)/(.......)(\d).(\w+).\w+.(\d\d\d\d)(\d\d)(\d\d).MT.fasta.gz!) {
        $aLabel  = $1;
        $sTag    = $2;
        $sNum    = $3;
        $prialt  = "mito";
        $date    = "$5-$6-$7";
    }

    else {
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

    $$data{"${prialt}${sNum}"}         = "${aLabel}${sNum}";
    $$data{"${prialt}${sNum}filesize"} = sprintf("%.3f GB", $filesize / 1024 / 1024 / 1024);
    $$data{"${prialt}${sNum}seq"}      = "https://s3.amazonaws.com/genomeark/$filename";

    ($$data{"${prialt}${sNum}n50ctg"},
     $$data{"${prialt}${sNum}n50scf"},
     $$data{"${prialt}${sNum}sizes"})  = generateAssemblySummaryHTML($prialt, $filename, $$data{"genome_size"});

    #  Update the assembly status based on the primary n50 and/or curation status.

    if ($prialt eq "pri") {
        if (($$data{"${prialt}${sNum}n50ctg"} < $goodCTG) ||
            ($$data{"${prialt}${sNum}n50scf"} < $goodSCF)) {
            $$data{"assembly_status"} = "<em style=\"color:red\">bad assembly</em>";
        }

        elsif ($aLabel !~ m/curated/) {
            $$data{"assembly_status"} = "<em style=\"color:orange\">preliminary assembly</em>";
        }

        else {
            $$data{"assembly_status"} = "<em style=\"color:green\">curated assembly</em>";
        }
    }
}



sub downloadAndSummarize ($$$) {
    my $name  = shift @_;
    my $size  = shift @_;
    my $file  = shift @_;
    my $bases = 0;

    if ((! -e "downloads/$name") &&
        (! -e "$name.summary")) {
        printf STDERR "FETCH file #%4d size %6.3f GB '%s'\n", $file, $size / 1024 / 1024 / 1024, $name;
        printf STDERR " -> downloads/$name\n";
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

    open(LS, "< genomeark.ls");
    while (<LS>) {
        chomp;

        my ($filedate, $filetime, $filesize, $filename) = split '\s+', $_;

        next   if ($filename !~ m!$name!);
        next   if ($filename =~ m!intermediate!);
        next   if ($filename =~ m!transcriptomic_data!);

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

        if    (($type eq "pacbio") && ($filename =~ m!genomic_data/pacbio/.*subreads.*bam$!)) {
            push @files, "$filesize\0$filename";
        }

        if    (($type eq "phase") && ($filename =~ m!genomic_data/phase/.*q.gz$!)) {
            push @files, "$filesize\0$filename";
        }
    }
    close(LS);

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

    return($scaling);
}



sub computeBionanoBases ($) {
    my $name    = shift @_;
    my @files;

    open(LS, "< genomeark.ls");
    while (<LS>) {
        chomp;

        my ($filedate, $filetime, $filesize, $filename) = split '\s+', $_;

        next   if ($filename !~ m!$name!);
        next   if ($filename =~ m!intermediate!);
        next   if ($filename =~ m!transcriptomic_data!);

        if    ($filename =~ m!genomic_data/bionano/.*bnx.gz$!) {
            push @files, "$filesize\0$filename";
        }
    }
    close(LS);

    #  Download all the cmap files.

    foreach my $f (@files) {
        my ($size, $name) = split '\0', $f;

        if ((! -e "downloads/$name") &&
            (! -e "$name.summary")) {
            printf STDERR "FETCH size %6.3f GB '%s'\n", $size / 1024 / 1024 / 1024, $name;
            system("aws --no-sign-request s3 cp s3://genomeark/$name downloads/$name")  if ($SKIP_RAW == 0);
        }
    }

    #  Parse each cmap to compute length of each contig.  Save.

if (0) {
    foreach my $f (@files) {
        my ($size, $name) = split '\0', $f;

        next  if (! -e "downloads/$name");
        next  if (  -e "$name.summary");

        if ($name =~ m/gz$/) {
            open(B, "gzip -dc downloads/$name |");
        } else {
            open(B, "< downloads/$name");
        }

        my %contigLen;

        while (<B>) {
            next  if (m/^#/);

            my @v = split '\s+', $_;

            $contigLen{$v[0]} = $v[1];
        }

        close(B);

        system("mkdir -p $name");
        system("rmdir    $name");

        open(B, "> $name.summary") or die "Failed to create '$name.summary': $!\n";
        print B "CMapId\tContigLength\n";

        foreach my $k (keys %contigLen) {
            print B "$k\t$contigLen{$k}\n";
        }

        close(B);
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

if (! -e "vgp-metadata") {
    print STDERR "FETCHING METADATA.\n";
    system("git clone git\@github.com:VGP/vgp-metadata.git");
}

if (! -e "genomeark.ls") {
    print STDERR "FETCHING AWS FILE LIST.\n";
    system("aws --no-sign-request s3 ls --recursive s3://genomeark/ > genomeark.ls");

    print STDERR "UPDATING METADATA.\n";
    system("cd vgp-metadata ; git fetch ; git merge");
}

if (! -e "genomeark.ls") {
    die "ERROR: no 'genomeark.ls' file list, can't update.\n";
}

my $lastUpdate = (stat("genomeark.ls"))[9];

my @speciesList = discover(@ARGV);
my %asmToShow   = loadAssemblyStatus();


#  Rebuild each species.md file.

foreach my $species (@speciesList) {
    my %seqFiles;
    my %seqBytes;
    my %meta;
    my %data;

    my $name;
    my $asm;

    loadMeta($species, \%meta);
    loadData($species, \%data);
    #printData($species, \%data);

    $name = makeName($meta{"species.name"});


    die "No meta{species.name} found?\n"  if ($name eq "");

    #
    #  Delete unused stuff.
    #

    delete $data{'status'};

    #
    #  Delete stuff we generate.
    #

    undef %data;

    #
    #  Copy metadata to the page.
    #

    $data{"name"}                = $meta{"species.name"};
    $data{"common_name"}         = $meta{"species.common_name"};
    $data{"taxon_id"}            = $meta{"species.taxon_id"};

    $data{"s3"}                  = "s3://genomeark/species/$name";

    $data{"genome_size"}         = $meta{"species.genome_size"};
    $data{"genome_size_display"} = prettifyBases($meta{"species.genome_size"});  #  Updated later, too.
    $data{"genome_size_method"}  = $meta{"species.genome_size_method"};

    $data{"image"}               = $meta{"species.image"};
    $data{"image_license"}       = $meta{"species.image_license"};

    $data{"data_status"}         = "<em style=\"color:red\">no data</em>";
    $data{"assembly_status"}     = "<em style=\"color:red\">no assembly</em>";

    #$data{"last_raw_data"}       = Do not set here; should only be present if data exists.
    $data{"last_updated"}        = $lastUpdate;

    #
    #  Find which assembly to use.
    #

    $data{"assembly"} = $asm    = $asmToShow{$name};

    print STDERR "\n";
    print STDERR "$name -- $asm\n";

    open(LS, "< genomeark.ls");
    while (<LS>) {
        chomp;

        my ($filedate, $filetime, $filesize, $filename) = split '\s+', $_;

        next   if ($filename !~ m!$name!);
        next   if ($filename =~ m!intermediate!);
        next   if ($filename =~ m!transcriptomic_data!);

        if ($filename =~ m!genomic_data!) {
            if ("$filedate $filetime" =~ m/(\d\d\d\d)-(\d\d)-(\d\d)\s+(\d\d):(\d\d):(\d\d)/) {
                my ($yr, $mo, $da, $hr, $mn, $sc) = ($1, $2, $3, $4, $5, $6);

                my $t = timelocal($sc, $mn, $hr, $da, $mo-1, $yr); 

                if ($data{"last_raw_data"} < $t) {     #  If this isn't set, Raw Data shows
                    $data{"last_raw_data"} = $t;       #  "No data.".
                }
            }

            processData($filesize, $filename, \%seqFiles, \%seqBytes);
            next;
        }

        #print "  '$filename'\n";

        if (($asm ne "") &&
            ($filename =~ m!$asm.*fasta!)) {
            #print "  '$filename' -- MATCH\n";
            processAssembly($filesize, $filename, \%data);
            next;
        }
    }
    close(LS);

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
    #printf "%14d gzip bytes in %3d Phase datasets.\n",             $seqBytes{"phase"} / 2,    $seqFiles{"phase"}       if ($seqFiles{"phase"} > 0);

    $seqFiles{"dovetail"} /= 2;
    $seqFiles{"illumina"} /= 2;
    $seqFiles{"phase"}    /= 2;

    foreach my $type (qw(10x arima bionano dovetail illumina pbscraps pbsubreads phase)) {
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
        $data{"data_pbsubreads_scale"} = estimateRawDataScaling($name, "pacbio");
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

    if (($seqBytes{"10x"} > 0)) {
        if ($data{"genome_size"} > 0) {
            $data{"data_10x_coverage"}        = sprintf("%.2fx", $seqBytes{"10x"} * $data{"data_10x_scale"} / $data{"genome_size"});
            $data{"data_10x_bases"}           = prettifyBases($seqBytes{"10x"} * $data{"data_10x_scale"});
        }
        $data10x++;
    }

    if (($seqBytes{"arima"} > 0)) {
        if ($data{"genome_size"} > 0) {
            $data{"data_arima_coverage"}      = sprintf("%.2fx", $seqBytes{"arima"} * $data{"data_arima_scale"} / $data{"genome_size"});
            $data{"data_arima_bases"}         = prettifyBases($seqBytes{"arima"} * $data{"data_arima_scale"});
        }
        $dataHIC++;
    }

    if (($seqBytes{"bionano"} > 0)) {
        if ($data{"genome_size"} > 0) {
            my $b = computeBionanoBases($name);
            $data{"data_bionano_coverage"}    = sprintf("%.2fx", $b / $data{"genome_size"});
            $data{"data_bionano_bases"}       = prettifyBases($b);
        }
        $dataBio++;
    }

    if (($seqBytes{"dovetail"} > 0)) {
        if ($data{"genome_size"} > 0) {
            $data{"data_dovetail_coverage"}   = sprintf("%.2fx", $seqBytes{"dovetail"} * $data{"data_dovetail_scale"} / $data{"genome_size"});
            $data{"data_dovetail_bases"}      = prettifyBases($seqBytes{"dovetail"} * $data{"data_dovetail_scale"});
        }
        $dataHIC++;
    }

    if (($seqBytes{"illumina"} > 0)) {
        if ($data{"genome_size"} > 0) {
            $data{"data_illumina_coverage"}   = sprintf("%.2fx", $seqBytes{"illumina"} * $data{"data_illumina_scale"} / $data{"genome_size"});
            $data{"data_illumina_bases"}      = prettifyBases($seqBytes{"illumina"} * $data{"data_illumina_scale"});
        }
    }

    if (($seqBytes{"pbsubreads"} > 0)) {
        if ($data{"genome_size"} > 0) {
            $data{"data_pbsubreads_coverage"} = sprintf("%.2fx", $seqBytes{"pbsubreads"} * $data{"data_pbsubreads_scale"} / $data{"genome_size"});
            $data{"data_pbsubreads_bases"}    = prettifyBases($seqBytes{"pbsubreads"} * $data{"data_pbsubreads_scale"});
        }
        $dataPac++;
    }

    if (($seqBytes{"phase"} > 0)) {
        if ($data{"genome_size"} > 0) {
            $data{"data_phase_coverage"}      = sprintf("%.2fx", $seqBytes{"phase"} * $data{"data_phase_scale"} / $data{"genome_size"});
            $data{"data_phase_bases"}         = prettifyBases($seqBytes{"phase"} * $data{"data_phase_scale"});
        }
        $dataHIC++;
    }

    if (($dataPac > 0) ||
        ($data10x > 0) ||
        ($dataHIC > 0) ||
        ($dataBio > 0)) {
        $data{"data_status"} = "<em style=\"color:orange\">some data</em>";
    }

    if (($dataPac > 0) &&
        ($data10x > 0) &&
        ($dataHIC > 0) &&
        ($dataBio > 0)) {
        $data{"data_status"} = "<em style=\"color:green\">all data</em>";
    }



    #printData($species, \%data);
    saveData($species, \%data);
}
