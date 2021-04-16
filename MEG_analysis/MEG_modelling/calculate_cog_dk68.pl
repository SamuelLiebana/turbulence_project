#!/usr/bin/perl

# calculate the centers of gravity for all 68 regions
# for 1:68
#   create binary mask from dk68 file provided by first argument
#   run fslstats -c and outputs to dk68_cog.txt

die "usage: calculate_cog_dk68.pl dk68\n" if @ARGV <1;
$dk68=$ARGV[0];

    $arg="rm dk68_cog.txt";
    print($arg."\n");
    system($arg);

for ($t=1; $t<69; $t++){

    $arg="fslmaths $dk68 -thr $t -uthr $t -bin $t";
    print($arg."\n");
    system($arg);

    $arg="fslstats $t.nii.gz -c >>dk68_cog.txt";
    print($arg."\n");
    system($arg);

    $arg="rm $t.nii.gz";
    print($arg."\n");
    system($arg);
    
} 

