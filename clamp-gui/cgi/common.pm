package common;
use strict; use warnings;
use CGI qw(:all *table *Tr *td *blockquote -nosticky);
use CGI::Carp qw(fatalsToBrowser);  autoEscape(undef);
use Exporter 'import';
our @EXPORT = qw($SITE $DESC ShowHead ShowFoot Fail $WADM $WNAM); 

$ENV{'PATH'} = '/bin:/usr/bin';
delete @ENV{'IFS', 'CDPATH', 'ENV', 'BASH_ENV'};

my %sites = ('cml', 'Philadelphia Variants',
	     'tp53','Tumor Protein p53',
	     'cmldev', 'CML Dev Server');
our $WADM = 'wesley.schaal@farmbio.uu.se';
our $WNAM = 'wes';

our ($SITE, $DESC); 
my ($root, $conf) = ('.', '');
my @dir = split(m#/#,$ENV{SCRIPT_NAME});
if ($dir[-2] eq 'reg') {
    $SITE  = $dir[-3];
    $root .= '.';
} else {
    $SITE = $dir[-2];
    $conf = 'reg/';
}
#our ($SITE) = $ENV{SCRIPT_NAME} =~ m#/([^/]+)/[^/]+$#;
#our $DESC = $sites{$SITE}; 
$DESC = $sites{$SITE}; 
$SITE = uc $SITE;
$DESC ||= "$SITE Project";

my @link = ('Admin', "${conf}config.cgi");

sub ShowHead {
    my @lnx = (@link,@_); 
    Fail("Header trouble: @lnx") unless $#lnx % 2;
    my ($v, $k);
    #$headed++;

    $\="\n";
    print header();
    print start_html(-title=>"$DESC",
		          author=>$WADM,
		          'bgcolor="#FFFFFF"');
    print '<center><font color="555555">[';
    print a({href=>$root}, $SITE);
    while($k = shift(@lnx)) {
	Fail("Header trouble: @lnx") unless $v = shift(@lnx);
	print ' | ' . a({href=>"$v"},"$k");
    }
    print ' ]</font>';
    print '<h3><img src="img/farmbio.gif"';
    print 'alt="<logo>"></h3></center>';
}


sub ShowFoot {
    my @lnx = (@link,@_); Fail("Footer trouble: @lnx") unless $#lnx % 2;
    my ($v, $k, @m, $m);
    @m = (localtime((stat($0))[9]))[3..5];
    $m = sprintf("%04d-%02d-%02d",$m[2]+1900,$m[1]+1,$m[0]);
    print br({clear=>'all'}),hr({noshade=>1});
    my $url = url(); # $url =~ s/http:/https:/ if https();
    $url =~ s|http://export.uppmax|https://export.uppmax|; # KLUDGE
    print font({size=>'-1'},"URL: $url" .
	       br . 'Copyright &copy; 2015 UAS/UU. All rights reserved.' .
	       br . "Last modified: $m by " .
	       a({href=>"mailto:$WADM"},$WNAM) .
	       #br . a({href=>'contacts', target=>'_new'},'Contacts') .
	       br .'Nav: [ ' .
	       a({href=>$root}, $SITE));
    while($k = shift(@lnx)) {
	Fail("Footer trouble: @lnx") unless $v = shift(@lnx);
	print font({size=>'-1'}," | " . a({href=>"$v"},"$k"));
    }
    print font({size=>'-1'},' ] ');
    print end_html;
}


sub Fail {
    #ShowHead(); # check for dup?
    my ($reason) = @_;
    print h2("Error: &nbsp; $reason.");
    print p('Please contact ' .
	        a({href=>"mailto:$WADM?subject=Error on " . url(-absolute=>1) .
		          '&body=Time: ' . localtime() . '%0AError: ' .
			     $reason . '%0A%0A'},$WNAM) .
	        ' if you think you have found an error.');
    print p('Return to ' . a({href=>referer()},'Originating Page') . "\?\n");
    print end_html;
    exit 0;
}

1;
