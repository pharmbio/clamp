#!/usr/bin/perl -T
#
# Name: Register samples/patients
# Auth: Wesley Schaal, wesley.schaal@farmbio.uu.se
# When: 2014-10-02

use lib ('../perl/lib/perl5', '..');
use v5.10; use strict; use warnings;
use CGI qw(:all *table *Tr *td *blockquote -nosticky);
use CGI::Carp qw(fatalsToBrowser);  autoEscape(undef);
use common; # local hacks
#require "./common.pl";
#$appl = 'Philadelphia Portal';
my $appl = lc $SITE; # common
my $desc = $DESC;    # common

#my $dbf = "../data/${appl}.sqlite";
#my $dbh = DBI->connect("dbi:SQLite:dbname=$dbf",'','');

ShowHead();
Chooser();
ShowFoot();


sub Chooser {
    print start_blockquote, h2("Register samples for $desc");

    # print p('Choose the task you wish to administrate:');
    print ul(li([
		 # a({href=>'register.cgi'}, 'Register a Sample') . br . br,
		 a({href=>'bulkreg.cgi'},  'Register Samples') . br . br,
		 a({href=>'loadfile.cgi'}, 'Upload Files') . br . br,
		 # a({href=>'editreg.cgi'}, 'Edit Registration') . br . br,
		 font({color=>'#C0C0C0'},'Edit Registration') . br . br,
		]));
    print end_form, end_blockquote;
}


sub Fret { 
    my $v = $_[0];
    return unless $v;
    print p(font({size=>'+1',color=>'green'},$v));
}
