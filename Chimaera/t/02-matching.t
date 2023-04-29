#!perl -T
use 5.006;
use strict;
use warnings;
use Test::More tests => 10;
use Test::Exception;


my $HAPLO_1 = "aaaaaaaCaaaaCaaaCaaaa";
my $HAPLO_2 = "aaaaaaaTaaaaTaaaTaaaa";

my $SHORTSQ = "aaaaaaaCaaaaCaaaTaaa";
my $LONGSEQ = "aaaaaaaCaaaaCaaaTaaaaa";

my $TESTccT = "aaaaaaaCaaaaTaaaTaaaa";
my $TESTcTT = "aaaaaaaCaaaaTaaaTaaaa";
my $TESTTcc = "aaaaaaaTaaaaCaaaCaaaa";
my $TESTTTc = "aaaaaaaTaaaaTaaaCaaaa";

my $TESTcTc = "aaaaaaaCaaaaTaaaCaaaa";
my $TESTTcT = "aaaaaaaTaaaaCaaaTaaaa";

my $CASEccT = "aaAAAaacaaaaTaaataAaa";

BEGIN {
    use_ok( 'Chimaera::Matcher' ) || print "Can't find the Chimaera::Matcher module!\n";
}


# Build a Chimaera::Matcher object
my $matcher = Chimaera::Matcher->new(

	'haplotype1' => $HAPLO_1,
	'haplotype2' => $HAPLO_2

);

ok( ! $matcher->possible_chimaera( $SHORTSQ ), "02 - NOT a possible chimaera - too short" );
ok( ! $matcher->possible_chimaera( $LONGSEQ ), "03 - NOT a possible chimaera - too long" );
ok( $matcher->possible_chimaera( $TESTccT ), "04 - Should report a possible chimaera" );
ok( $matcher->possible_chimaera( $TESTcTT ), "05 - Should report a possible chimaera" );
ok( $matcher->possible_chimaera( $TESTTcc ), "06 - Should report a possible chimaera" );
ok( $matcher->possible_chimaera( $TESTTTc ), "07 - Should report a possible chimaera" );

ok( ! $matcher->possible_chimaera( $TESTcTc ), "08 - NOT a possible chimaera - multiple crossovers" );
ok( ! $matcher->possible_chimaera( $TESTTcT ), "09 - NOT a possible chimaera - multiple crossovers" );

ok( $matcher->possible_chimaera( $CASEccT ), "10 - Should report a possible chimaera" );
