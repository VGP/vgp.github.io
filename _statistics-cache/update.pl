#!/usr/bin/perl

use strict;





#
#  Discover species.
#

sub discover {
    my @speciesList;

    open(LS, "ls ../_genomeark/*.md |");
    while (<LS>) {
        chomp;

        if (m!genomeark/(\w+).md$!) {
            #print STDERR "Species: '$1'\n";
            push@speciesList, $1;
        }
    }
    close(LS);

    return(@speciesList);
}


#
#  Load the existing .md page for a species.
#

sub loadData ($$$) {
    my $species = shift @_;
    my $keys    = shift @_;
    my $data    = shift @_;

    my  $key = undef;

    undef @$keys;
    undef %$data;

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
            push @$keys, $key;
            $$data{$key} = "|\n";
        }

        #  Add data to multi-line parameters.
        elsif (defined($key) && (m!^\s\s\S!)) {
            $$data{$key} .= "$_\n";
        }

        #  Match any key-value pairs.
        elsif (m!(\w+):\s*(.*)$!) {
            $key      = undef;
            push @$keys, $1;
            $$data{$1} = $2;
        }

        #  Match any extra stuff in the file.
        else {
            $$data{":"} .= "$_\n";
        }
    }
    close(MD);
}



sub saveData ($$$) {
    my $species = shift @_;
    my $keys    = shift @_;
    my $data    = shift @_;

    #  Scan the keys used in data, add any new ones to the end of
    #  the key list.

    foreach my $k (keys %$data) {
        foreach my $s (@$keys) {
            if ($k eq $s) {
                undef $k;
                last;
            }
        }

        if (defined($k)) {
            push @$keys, $k;
        }
    }

    open(MD, "> ../_genomeark/$species.md");

    @$keys = sort @$keys;

    print MD "---\n";
    foreach my $key (@$keys) {
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
    my $keys    = shift @_;
    my $data    = shift @_;

    #foreach my $key (@$keys) {
    #    print STDERR "KEY='$key'\n";
    #}

    foreach my $key (@$keys) {
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

    #  Shows up to "### bp"
    if      ($size < 1000) {
        return("$size  bp");
    }

    #  Shows up to "##.## Kbp"
    elsif ($size < 100000) {
        return(sprintf("%.2f Kbp", $size / 1000));
    }

    #  Shows up to "##.## Mbp"
    elsif ($size < 100000000) {
        return(sprintf("%.2f Mbp", $size / 1000000));
    }

    #  Everything else.  I suspect this code path will never be tested.
    else {
        return(sprintf("%.2f Gbp", $size / 1000000000));
    }
}



sub generateSizesHTML ($$) {
    my $filename   = shift @_;
    my $genomeSize = shift @_;
    my $n50table;
    my $n50;

    $filename =~ s/.gz//;

    if ((! -e "$filename.gz") &&
        (! -e "$filename.summary")) {
        print STDERR "Fetch s3://genomeark/$filename.gz\n";

        system("mkdir -p $filename");
        system("rmdir    $filename");
        system("aws --no-sign-request s3 cp s3://genomeark/$filename.gz $filename.gz");
    }

    if (! -e "$filename.summary") {
        print STDERR "SUMMARIZING '$filename' with genome_size $genomeSize\n";
        system("sequence summarize -size $genomeSize $filename.gz > $filename.summary");
    }

    if (! -e "$filename.summary") {
        print STDERR "FAILED TO FIND SIZES '$filename'\n";
        return(undef);
    }


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
    $n50table .= "  <th>NG</th>\n";
    $n50table .= "  <th>LG</th>\n";
    $n50table .= "  <th>Length</th>\n";
    $n50table .= "  <th>Cumulative</th>\n";
    $n50table .= "  </tr>\n";
    $n50table .= "  </thead>\n";
    $n50table .= "  <tbody>\n";

    #  Strip out the N50 chart

    open(SU, "< $filename.summary");
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

        if ($a == 0) {   #  If zero, the column was empty, or not a number.
            next;
        }

        if ($a =~ m!x!) {
            $a =~ s/^0+//;      #  Strip off leading zeros from "000.851x"
            $a =~ s/^\./0./;    #  But then add back one to get "0.851x"
            $c = prettifyBases($c);

            $n50table .= "  </tbody>\n";
            $n50table .= "  <tfoot>\n";
            $n50table .= "  <tr><th>$a</th><th>$b</th><th></th><th>$c</th></tr>\n";
            $n50table .= "  </tfoot>\n";
            last;
        }

        else {
            $a =~ s/^0+//;
            $b = prettifyBases($b);
            $d = prettifyBases($d);

            $n50table .="  <tr><td>$a</td><td>$c</td><td>$b</td><td>$d</td>\n";
        }
    }
    close(SU);

    $n50table .= "  </table>\n";
}




sub processData ($$$$) {
    my $filesize = shift @_;
    my $filename = shift @_;
    my $seqFiles = shift @_;
    my $seqBytes = shift @_;

    #  Parse out species and individual.

    #my $sName;
    #my $sTag;
    #my $sNum;
    #
    #if ($filename =~ m!^species/(.*_.*)/(.......)([0-9])/!) {
    #    $sName = $1;
    #    $sTag  = $2;
    #    $sNum  = $3;
    #
    #    $sTags{$sTag}++;
    #} else {
    #    print STDERR "Failed to parse species/individual from '$_'\n";
    #}

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
    my $filename = shift @_;
    my $keys     = shift @_;
    my $data     = shift @_;

    print "PROCESS $_\n";

    my ($aLabel, $sTag, $sNum, $prialt, $date) = undef;

    if    (m!assembly_(.+)/(.......)(\d).(\w+).\w+.(\d\d\d\d)(\d\d)(\d\d).fasta.gz!) {
        $aLabel  = $1;
        $sTag    = $2;
        $sNum    = $3;
        $prialt  = $4;
        $date    = "$5-$6-$7";
    }

    elsif (m!assembly_(,+)/(.......)(\d).(\w+).\w+.(\d\d\d\d)(\d\d)(\d\d).MT.fasta.gz!) {
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

    #  Enable this line to print the genomeark.ls like for the any assembly we're going to process.
    print "$_\n";

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

        #$$data{"tag${sNum}"}           = undef;

        #$$data{"pri${sNum}seq"}        = undef;
        #$$data{"pri${sNum}sizes"}      = undef;

        #$$data{"alt${sNum}seq"}        = undef;
        #$$data{"alt${sNum}sizes"}      = undef;

        #$$data{"mito${sNum}seq"}       = undef;
        #$$data{"mito${sNum}sizes"}     = undef;
    }

    #  Update everything for this file.

    #print "  add prialt $prialt sNum ${sNum} date $date for $filename\n";

    $$data{"${prialt}${sNum}"}         = "${aLabel}${sNum}";
    $$data{"${prialt}${sNum}filesize"} = sprintf("%.3f GB", (-s $filename) / 1024 / 1024 / 1024);
    $$data{"${prialt}${sNum}seq"}      = "https://s3.amazonaws.com/genomeark/$filename";
    $$data{"${prialt}${sNum}sizes"}    = generateSizesHTML($filename, $$data{"genome_size"});
}



sub estimatePacBioScaling ($) {
    my $name    = shift @_;
    my @files;

    open(LS, "< genomeark.ls");
    while (<LS>) {
        chomp;

        my ($filedate, $filetime, $filesize, $filename) = split '\s+', $_;

        next   if ($filename !~ m!$name!);
        next   if ($filename =~ m!intermediate!);
        next   if ($filename =~ m!transcriptomic_data!);

        if    ($filename =~ m!genomic_data/pacbio/.*subreads.*bam$!) {
            push @files, "$filesize\0$filename";
        }
    }
    close(LS);

    @files = sort { $a <=> $b } @files;

    my $nFiles = scalar(@files);

    my $f1 = int(1 * $nFiles / 4);
    my $f2 = int(2 * $nFiles / 4);
    my $f3 = int(3 * $nFiles / 4);

    my ($size1, $name1, $bases1, $seqs1) = split '\0', @files[$f1];
    my ($size2, $name2, $bases2, $seqs2) = split '\0', @files[$f2];
    my ($size3, $name3, $bases3, $seqs3) = split '\0', @files[$f3];

    if ((! -e "$name1") &&
        (! -e "$name1.summary")) {
        printf STDERR "FETCH file #%4d size %6.3f GB '%s'\n", $f1, $size1 / 1024 / 1024 / 1024, $name1;
        system("aws --no-sign-request s3 cp s3://genomeark/$name1 $name1");
    }
    if ((  -e "$name1") &&
        (! -e "$name1.summary")) {
        printf STDERR "EXTRACT bases\n";
        system("samtools bam2fq $name1 > $name1.fastq");
        printf STDERR "SUMMARIZE $name1.summary\n";
        system("sequence summarize $name1.fastq > $name1.summary");
    }
    if (  -e "$name1.summary") {
        open(ST, "< $name1.summary");
        while (<ST>) {
            if (m/^G=(\d+)\s/) {
                $bases1 = $1;
                last;
            }
        }
        close(ST);
    }

    if ((! -e "$name2") &&
        (! -e "$name2.summary")) {
        printf STDERR "FETCH file #%4d size %6.3f GB '%s'\n", $f2, $size2 / 1024 / 1024 / 1024, $name2;
        system("aws --no-sign-request s3 cp s3://genomeark/$name2 $name2");
    }
    if ((  -e "$name2") &&
        (! -e "$name2.summary")) {
        printf STDERR "EXTRACT bases\n";
        system("samtools bam2fq $name2 > $name2.fastq");
        printf STDERR "SUMMARIZE $name2.summary\n";
        system("sequence summarize $name2.fastq > $name2.summary");
    }
    if (  -e "$name2.summary") {
        open(ST, "< $name2.summary");
        while (<ST>) {
            if (m/^G=(\d+)\s/) {
                $bases2 = $1;
                last;
            }
        }
        close(ST);
    }

    if ((! -e "$name3") &&
        (! -e "$name3.summary")) {
        printf STDERR "FETCH file #%4d size %6.3f GB '%s'\n", $f3, $size3 / 1024 / 1024 / 1024, $name3;
        system("aws --no-sign-request s3 cp s3://genomeark/$name3 $name3");
    }
    if ((  -e "$name3") &&
        (! -e "$name3.summary")) {
        printf STDERR "EXTRACT bases\n";
        system("samtools bam2fq $name3 > $name3.fastq");
        printf STDERR "SUMMARIZE $name3.summary\n";
        system("sequence summarize $name3.fastq > $name3.summary");
    }
    if (  -e "$name3.summary") {
        open(ST, "< $name3.summary");
        while (<ST>) {
            if (m/^G=(\d+)\s/) {
                $bases3 = $1;
                last;
            }
        }
        close(ST);
    }

    if (($bases1 == 0) ||
        ($bases2 == 0) ||
        ($bases3 == 0)) {
        die "FAILED TO ESTIMATE SIZES.  $bases1 $bases2 $bases3\n";
    }

    my $scaling = ($bases1 + $bases2 + $bases3) / ($size1 + $size2 + $size3);

    return($scaling);
}



#
#  Main
#

if (! -e "genomeark.ls") {
    system("aws --no-sign-request s3 ls --recursive s3://genomeark/ > genomeark.ls");
}

my @speciesList = discover();

if (scalar(@ARGV > 0)) {
    @speciesList = @ARGV;
}

#  Check that everything is sane.

my $notSane = 0;

foreach my $species (@speciesList) {
    my @keys;
    my %data;

    loadData($species, \@keys, \%data);

    if ($data{"genome_size"} == 0) {
        print STDERR "$species has no genome_size property.\n";
        $notSane++;
    }
}

die "Problems!\n"  if ($notSane);

#  Rebuild each species.md file.

foreach my $species (@speciesList) {
    my %seqFiles;
    my %seqBytes;
    my @keys;
    my %data;

    loadData($species, \@keys, \%data);
    #printData($species, \@keys, \%data);

    my $name = $data{"name"};
    $name =~ s/\s+/_/g;

    my $asm  = $data{"assembly"};
    $asm = "no-assembly-to-show"   if ($asm eq "");

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
            processData($filesize, $filename, \%seqFiles, \%seqBytes);
            next;
        }

        print "  '$filename'\n";

        if ($filename =~ m!$asm!) {
            print "  '$filename' -- MATCH\n";
            processAssembly($filename, \@keys, \%data);
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
            $data{"data_${type}_bytes"}    = sprintf("%.3f GB", $seqBytes{$type} / 1024 / 1024 / 1024);
            $data{"data_${type}_coverage"} = "<em style=\"color:red\">unknown</em>";
            $data{"data_${type}_bases"}    = "<em style=\"color:red\">unknown</em>";
            $data{"data_${type}_files"}    = $seqFiles{$type};
        }
    }

    if ($data{"data_pbsubreads_scale"} == 0) {
        $data{"data_pbsubreads_scale"} = estimatePacBioScaling($name);
    }

    die "Invalid genome_size?\n"  if ($data{"genome_size"} == 0);  #  Checked before we start this, so should never occur.

    $data{"data_pbsubreads_coverage"} = sprintf("%.2fx", $seqBytes{"pbsubreads"} * $data{"data_pbsubreads_scale"} / $data{"genome_size"});
    $data{"data_pbsubreads_bases"}    = prettifyBases($seqBytes{"pbsubreads"} * $data{"data_pbsubreads_scale"});

    #$seqFiles += $seqFiles{"pbsubreads"};
    #$seqBytes += $seqBytes{"pbsubreads"};

    #printData($species, \@keys, \%data);
    saveData($species, \@keys, \%data);
}
