#!/usr/bin/perl

$PRECISION = 6;

my $meta = 1;
my $n_meta = 0;
my $n_data = 0;

# Reads from stdin
# Expected format, two lines per record:
# iteration algo, layer, is_fwd
# set of 12 float parameters (or any number, as expected by fill_hit_selection_windows_params())

while ($_ = <>)
{
    chomp $_;
    my @arr = split(/\s+/, $_);
    if ($meta) {
        push @meta, @arr;
        ++$n_meta;
    } else {
        map { $_ = sprintf "%.${PRECISION}g", $_;
              $_ =~ s/(.*\..+)0+(e[-+]\d\d)?$/$1$2/;
              $_ = "0.0" if $_ eq "0";
              $_ } @arr;
        push @data, @arr;
        ++$n_data;
    }
    $meta = 1 - $meta;
}

if ($n_meta != $n_data) {
    die "Different number of meta / data lines."
} else {
    print "Read $n_meta entries for meta and data ... creating .h file.\n";
}

open OH, ">CMS-phase1-HitSelectionWindows.h" or die "Can not open output file";

print OH "std::vector<int> meta = { ", join(", ", @meta), "};\n";
print OH "std::vector<float> data = { ", join(", ", @data), "};\n";

close OH;
