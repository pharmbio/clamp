#!/usr/bin/perl
#
# Name: Register samples/patients
# Auth: Wesley Schaal, wesley.schaal@farmbio.uu.se
# When: 2015-01-14

use lib ('../perl/lib/perl5', '..');
use DBI; use strict; use warnings;
use CGI qw(:all *table *Tr *td *blockquote -nosticky);
use CGI::Carp qw(fatalsToBrowser);  autoEscape(undef);
use Time::Piece; use Time::Seconds;
use common; # local hacks
#require "./common.pl";
#$appl = 'Philadelphia Portal';
my $appl = lc $SITE; # common
my $desc = $DESC;    # common
my $data = '../uploads';

my $dbf = "../../data/${appl}.sqlite";
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbf",'','');

my ($D, $submit) = (0, param('submit'));
ShowHead(qw(Register bulkreg.cgi Upload loadfile.cgi));
ShowForm();
ShowFoot(qw(Register register.cgi));


sub ShowForm {
    my $topn = 3;
    print start_blockquote;
    print h2('Upload data files for registered samples');
    print p('File names determine the "Run ID" ' .
	    'they will be associated with. Use caution. ' .br. 'All files ' .
	    'are technically optional but it\'s rather boring without "Results".');

    print start_multipart_form, br;
    print start_table({border=>1,align=>'center',width=>'90%',
		       cellspacing=>0,cellpadding=>'6'});
    print Tr(td(['Results file (.txt)' . br .
		 filefield(-name=>'txt'),
		 'Sequences (.fastq.gz)' . br .
		 filefield(-name=>'seq'),
		 'Coverage (QC.pdf)' . br .
		 filefield(-name=>'pdf'),
		]));
    print Tr(td(['Clonal Distribution (.txt)' . br .
		 filefield(-name=>'dtab'),
		 'Clonal Distribution (.pdf)' . br .
		 filefield(-name=>'dist'),
		 'Log file (log.txt)' . br .
		 filefield(-name=>'log'),
		]));
#    print Tr(td(['Data files (txt, pdf or fastq.gz)' . br .
#		 filefield(-name=>'data')]));
    print end_table;
    print table({align=>'center',width=>'50%',cellspacing=>'10'},
		Tr(th([submit(-name=>'submit',value=>'Submit'),
		       defaults('Reset')])));
    print end_form;

    # print hr({width=>'75%'}), p('Recent registrations:');
    # my $sql = "select samid,runid,sdate from samples order by regid desc limit $topn";
    # my $sth = $dbh->prepare($sql); $sth->execute();
 
    # print start_blockquote, start_table({cellpadding=>5,border=>1});
    # print Tr(th(['Sample ID', 'Run ID', 'Date']));
    # while(my @row = $sth->fetchrow_array) {
    # 	print Tr(td([@row]));
    # }
    # print end_table, end_blockquote, end_blockquote;

    NewEntry() if $submit;
}


sub NewEntry {
    my $ruser = remote_user();
    my $rhost = remote_host();
#    my $sth = $dbh->prepare('insert into files (regid,versn,ftype,ruser,rhost) '.
#     		    'values (?,?,?,?,?)');
    my $sth = $dbh->prepare('insert into files (regid,versn,orig,ftype,ruser,rhost) values (?,?,?,?,?,?)');
    my $stc = $dbh->prepare("select regid from samples where runid = ?");
    my $stv = $dbh->prepare("select max(versn) from files where regid = ? and ftype = ?");

    # my @typ = qw(txt seq pdf dist);
    # my %typ; @typ{@typ} = qw(txt fastq.gz pdf txt);
    my @typ = qw(txt seq pdf dist dtab log);
    my @suf = qw(mutations_final.txt .fastq.gz QC.pdf clonal_distribution.pdf clonal_distribution.txt log.txt);
    my %typ; @typ{@typ} = @suf;
    
    print hr, start_blockquote;
    print start_table({cellpadding=>5,border=>1});
    print Tr(th([qw(Type Original Parsed RegID Version Status)]));

    for my $typ (@typ) {
     	my $ori = param($typ);
	print start_Tr, td([$typ, $ori]);
	unless ($ori) { 
	    print td(['&nbsp;','&nbsp;','&nbsp;','No file ... skipping']), end_Tr;
	    next;
	}

	my ($vic) = $ori =~ /([a-z_]+_\d+_\d+)/i;
	print td($vic);
	unless ($vic) {
	    print td(['&nbsp;','&nbsp;',"Can't get a Run ID from filename"]), end_Tr;
	    next;
	}
	# check suffix?

	my $regid = $dbh->selectrow_array($stc,undef,($vic));
	unless($regid) { 
	    print td(['&nbsp;','&nbsp;','No registered sample found']), end_Tr;
	    next;
	}
	# my $versn = $dbh->selectrow_array($stv,undef,($regid,$typ)) + 1;
	my $versn = $dbh->selectrow_array($stv,undef,($regid,$typ));

	my ($odir, $mkdir) = ("$data/$vic", '');
	unless(-e $odir) { mkdir("$odir") || ($mkdir = $!); }
	if ($mkdir) { 
	    print td(['&nbsp;','&nbsp;',"Error: can't create folder '$vic': $mkdir"]), end_Tr;
	    last;
	}

    	my $fil = upload($typ);
    	next unless defined $fil;
    	my $in = $fil->handle;
    	my $fn = sprintf "%s/%s", $odir, $ori;

	if($versn) {
	    $_='';
	    unless(rename($fn, sprintf("%s.%02d", $fn, $versn))) {
		print td(['&nbsp;','&nbsp;',"Error: can't archive prev '$vic': $_"]), end_Tr;
		last;
	    }
	}
	$versn++;  print p("OUT: $fn") if $D;

    	open (my $ut,'>',$fn);
    	local undef $\;
    	while ($in->read(my $buffer, 1024)) { print $ut $buffer; }
    	close $in; close $ut;
    	$sth->execute($regid, $versn, $ori, $typ, $ruser, $rhost);
	print td([$regid,$versn,'File uploaded']), end_Tr;
	ParseRes($fn, $regid, 0) if $typ eq 'txt';
	ParseLog($fn, $regid, 0) if $typ eq 'log';
    }
    print end_table;

    print p(font({color=>'darkgreen',size=>'+1'},'Unless errors were reported, '.
		 'your files should now be in the database and availble for viewing. '.
		 br. 'The system is ready for you to upload additional files.'));
    print end_blockquote;
    Delete_all();
    0;
}


sub ParseLog {
    my($fil, $regid, $bulk) = ($_[0], $_[1], $_[2]);
    my($log, $date, $type, $ref) = ('','','','');

    die("Bizzare regid '$regid'") unless $regid =~ /^\d+$/ && $regid > 0;
    print "\tParse: $regid: $fil\n" if $bulk;
    Fret("Empty log file '$fil'"),return(0) unless -s $fil;

    open(my $fh, '<', $fil) or Fret("Can't open log file '$fil'"),return(0);
    { local $/ = undef; $log = <$fh>; } 


    unless ($log =~ /CLAMP analysis complete/) {
	Fret('Incomplete Analysis'); 
    }
    if ($log =~ /^\w{3}\s+([^\-]+)\s+/) {
	$date = Time::Piece->strptime($1,"%b %d %T %Y");
	$date = $date->datetime;
    }
    if ($log =~ /Reference to be used for further analysis:\s\'([\w-]+)\'/) {
	$ref = $1; $ref =~ s/^\s*|\s*$//;
    } else {
	Fret("Unable to find ref for '$fil'");
	return(0);
    }
    if ($log =~ /Analysis type:\s*([\w\s-]+)$/m) {
	$type = $1; $type =~ s/^\s*|\s*$//;
    }
    print "\tref: $ref\t$type\t$date\n";
    # above regex make values safe for DB
    $dbh->do("insert into logs values ($regid, '$date', '$type', '$ref')");
}


sub ParseRes {
    #
    # Consider changing input parsing by header rather than column number
    # If so, warn or fail on added columns? Missing but unnecessary columns?
    #

    my $cols = "mutation,sequence,wt_reads_fwd,mut_reads_fwd,other_reads_fwd," .
	"freq_fwd,wt_reads_rev,mut_reads_rev,other_reads_rev,freq_rev," .
	"wt_reads,mut_reads,other_reads,freq,detection,routine";
    my $ncol = $cols =~ y/,/,/; $ncol++;

    my($fil, $regid, $bulk, $n) = ($_[0], $_[1], $_[2], 0); 
    die("Bizzare regid '$regid'") unless $regid =~ /^\d+$/ && $regid > 0;
    # print p("regid= $regid" .br. "fil = $fil") if $D;
    print "\tParse: $regid: $fil\n" if $bulk;

    Fret("Empty results file '$fil'"),return unless -s $fil;
    open(my $fh, '<', $fil) or Fret("Can't open results file '$fil'"),return(1);
    $_ = <$fh>; chomp; my $fld = split; s/^\s*|\s*$//g; s/\s+/,/g; my $head = $_;

    # print p("fld= $fld" .br. $head) if $D;
    Fret("$fld fields in results file ($ncol expected)"),return unless $fld == $ncol;
    Fret("Results header mismatch '$_' v '$head'"),return unless $head eq $cols;

    unless($bulk) {
	$dbh->do("delete from results where regid = $regid") or die($dbh->errstr());
    }

    my $sth = $dbh->prepare("insert into results (regid,sortby,$cols) " .
    			    'values (?,?' . ',?' x $ncol . ')');
    while(<$fh>) {
    	chomp; $n++;
    	$sth->execute($regid,$n,split) or Fret("DB error ($regid): $DBI::errstr");
    	print "$regid,$n,",split,"\n" if $D;
    }
    #print p("res = $n") if $D;
}


sub Fret {
    my $v = $_[0];
    return unless $v;
    print p(font({color=>'red'},$v));
}
