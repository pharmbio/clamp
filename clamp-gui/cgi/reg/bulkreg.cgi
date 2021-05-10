#!/usr/bin/perl -T
#
# Name: Register multiple samples/patients
# Auth: Wesley Schaal, wesley.schaal@farmbio.uu.se
# When: 2014-10-05

use lib ('../perl/lib/perl5','..');
use DBI; use strict; use warnings;
use CGI qw(:all *table *Tr *td *blockquote -nosticky);
use CGI::Carp qw(fatalsToBrowser);  autoEscape(undef);
use Time::Piece; use Time::Seconds;
use common; # local hacks
#require "./common.pl";
#my $WADM = 'wesley.schaal@farmbio.uu.se';
my $appl = lc $SITE; # common
my $desc = $DESC;    # common

#my $appl = 'Philadelphia Variants';

#my $dbf = '../../data/cml.sqlite';
my $dbf = "../../data/${appl}.sqlite";
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbf",'','');

# my $db   = $appl;
# # my $grp  = server_name() eq 'pele.farmbio.uu.se' ? "${appl}loc" : "${appl}bmc";
# my $grp  = "${appl}loc";
# my $dbh  = DBI->connect("DBI:mysql:$db"
# 	   . ";mysql_read_default_file=../data/.my.cnf"
# 	   . ";mysql_read_default_group=$grp") ||
#            die "Can't connect to $db: $DBI::errstr";

my @srcid = @{ $dbh->selectcol_arrayref('select land || "-" || site from sources order by srcid') };
my %srcid;  @srcid{@srcid} = (1 .. @srcid);

my @primr = @{ $dbh->selectcol_arrayref('select name from primers order by prmid') };
my %primr;  @primr{@primr} = (1 .. @primr);

my ($D, $status, $submit) = (0, 0, param('submit'));
ShowHead(qw(Register bulkreg.cgi));
if    ($submit eq 'Check') { CheckDat(); }
elsif ($submit eq 'Save')  { StoreDat(); }
else                       { ShowForm(); }
ShowFoot(qw(Register bulkreg.cgi Upload loadfile.cgi));


sub ShowForm {
    print h2("Register samples for $desc"),br;
    print start_blockquote;
    print p(b('Instructions:') .br.
	    'Enter tab- or space-separated text for ' . b('Run ID') . ', ' . 
	    b('Sample ID') . ', ' . b('Name of Primer') . ', ' . b('Site ID') .
	    ' and ' . b('Isolation Date') . '.' .br. 'These can be copy-pasted ' .
	    'from a spreadsheet or text file. Don\'t worry about exact spacing.' . br . 
	    'Dates must be in "yyyy-mm-dd" format. Missing dates mean "today". ' .
	    br. 'Quotes and other weird characters will be skipped. ' .br.
	    'Spaces within IDs will probably confuse the interpreter.' .br.br.
	    'After clicking the "Check" button below, you will get a page showing how the registration ' .br.
	    'was interpreted so you can confirm or reject before anything is saved.');
    print end_blockquote;

    print start_form;
    print start_table({width=>'120',cellspacing=>0,cellpadding=>0});
#    print start_table({align=>'center',width=>'99',cellspacing=>0,cellpadding=>0});
#    print Tr(td(i('Example:') . br . 
#		'pb_001_1 &nbsp; &nbsp; R1234 &nbsp; &nbsp; &nbsp; 2015-01-01'));
    print Tr(td(['&nbsp','Run Id', 'Sample Id', 'Primer', 'Site #', 'Date']));
    print Tr(td({colspan=>7}, textarea(-name=>'reg',rows=>13,cols=>110)));
    print Tr(td([i('Ex:'),qw(pb_001_1 R1234 ABL1_E1 UU 2015-01-01)]));
    print end_table,br;
    print table({width=>'50%',cellspacing=>'10'},
		Tr(th([submit(-name=>'submit',value=>'Check'),
		       defaults('Reset')])));
    print end_form;
}


sub CheckDat {
    print start_blockquote;
    print h2("Check registration data"),br;
    print p('The table below shows how the registration data will be interprepreted.');
    print ul(li('Samples already in the database are shown here with a grey ' .
		'background. These will be skipped if you try to save.'),
	     li('Uninterpretable or missing fields will be marked with a red ' .
		'background. These will be skipped if you try to save.'),
	     li('Suspicious fields will be marked with a green background. ' .
		'These are probably wrong but can be saved if you know they are correct.'),
	     li('Interpretation is rather strict but the details can be ' . 
		a({href=>"mailto:$WADM?subject=$SITE: reg format"},'negotiated') . '.')); 
    print p('Go back to the previous page to fix errors (if any) ' .
	    'or click ' . b('Save') . ' to complete the registration. ' .br.
	    'Note: Rows with grey or red cells will be silently skipped.');

    my (%run, %sam); my $t = localtime;
    my @color = ("#FFFFFF","#FFC0C0","#C0C0C0","#C9FFC9","C0C0F0");
    #my @primr = @{ $dbh->selectcol_arrayref('select name from primers order by prmid') };


    print start_table({align=>'center',width=>'95%',border=>1,cellspacing=>0,cellpadding=>6});
    print Tr(th(['Run ID', 'Sample ID', 'Primer', 'Source ID', 'Sample Date', 'Notes']));
    for my $row (split(/[\r\n]+/, param('reg'))) {
#    	my ($runid, $samid, $sdate) = split(' ',$row); 
    	my ($runid, $samid, $primr, $srcid, $sdate) = split(' ',$row); 
	my (@err, $not);
	$sdate ||= $t->ymd; # accept empty dates

	#my (@note, $note);

	($runid, $err[0]) = chk_runid($runid);
	($samid, $err[1]) = chk_samid($samid);
	($primr, $err[2]) = chk_primr($primr);
	($srcid, $err[3]) = chk_srcid($srcid);
	($sdate, $err[4]) = chk_sdate($sdate);

	$err[0] = 2 if $run{$runid}; $run{$runid}++;
	$err[1] = 2 if $sam{$samid}; $sam{$samid}++;

	($not, @err) = feedback(@err); $not ||= 'OK';

    	print start_Tr;
	print td({bgcolor=>$color[$err[0]]}, $runid);
	print td({bgcolor=>$color[$err[1]]}, $samid);
	print td({bgcolor=>$color[$err[2]]}, $primr);
	print td({bgcolor=>$color[$err[3]]}, $srcid);
	print td({bgcolor=>$color[$err[4]]}, $sdate);
	print td($not);
    	print end_Tr;
    }
    print end_table, br;
    print start_form;
    print table({align=>'center',width=>'95%',cellspacing=>'10'},
		Tr(th([submit(-name=>'submit',value=>'Save'),
		       defaults('Reset')])));
    print hidden('reg',param('reg'));
    print end_form;
    print end_blockquote;
}

sub StoreDat {
    my $t = localtime;
    my ($nrec, $ruser, $rhost, %run) = (0, remote_user(), remote_host());
    my $sav = $dbh->prepare('insert into samples (runid,samid,primr,srcid,'.
			    'sdate,ruser,rhost) values (?,?,?,?,?,?,?)') or die;

    # print p($dbh->selectrow_array('select count(*) from samples'));

    print h3('The following samples have been registered:');
    print start_table({border=>1,cellspacing=>0,cellpadding=>6});
#   print Tr(th([qw(Run Sample Date)]));
    print Tr(th(['Run ID', 'Sample ID', 'Primer', 'Source ID', 'Sample Date']));
    for my $row (split(/[\r\n]+/, param('reg'))) {
#    	my ($runid, $samid, $sdate) = split(' ',$row); 
    	my ($runid, $samid, $primr, $srcid, $sdate) = split(' ',$row); 
	$sdate ||= $t->ymd;
	my $err;
	
	($runid, $err) = chk_runid($runid);
	next if $err == 1 || $err == 2;
	next if $run{$runid}; $run{$runid}++;

	($samid, $err) = chk_samid($samid);
	next if $err == 1;

	($primr, $err) = chk_primr($primr);
	next if $err == 1;

	($srcid, $err) = chk_srcid($srcid);
	next if $err == 1;

	($sdate, $err) = chk_sdate($sdate);
	next if $err == 1;

#	$sav->execute($runid, $samid, $sdate, $ruser, $rhost) or die($dbh->errstr());
        $sav->execute($runid, $samid, $primr, $srcid, $sdate, $ruser, $rhost)
	    or die($dbh->errstr());
#    	print Tr(td([$runid,$samid,$sdate])); $nrec++;
    	print Tr(td([$runid, $samid, $primr, $srcid, $sdate])); $nrec++;
    }
    print end_table;
    print p("$nrec records stored");
    Delete_all();
}

sub chk_runid {
    my ($vic) = $_[0]; $vic =~ s/^\s*|\s$//g;
    my ($tmp) = $vic  =~ /([\w_-]+)/;
    return($vic, 1) unless $tmp;
    return($tmp, 1) unless $vic eq $tmp;
    return($tmp, 2) if $dbh->selectrow_array("select 1 from samples where runid= '$tmp'");
#   return($tmp, 3) unless $tmp =~ /^pb[\w_]+_\d+$/;
    return($tmp, 3) unless $tmp =~ /^\w+_\d+_\d+$/;
    return($tmp, 0);
}


sub chk_samid {
    my ($vic) = $_[0]; $vic =~ s/^\s*|\s$//g;
    my ($let,$tmp) = $vic  =~ /([A-Z]?)([\w_-]+)/;
    return($vic, 1) unless $tmp;
    return($vic, 1) unless $tmp < 1e8;
    $tmp = "${let}$tmp";
    return($tmp, 3) unless $let eq 'R';
    return($tmp, 3) unless $vic eq $tmp;
    return($tmp, 2) if $dbh->selectrow_array("select 1 from samples where samid= '$tmp'");
    return($tmp, 0);
}

sub chk_primr {
    my ($vic) = $_[0]; $vic =~ s/^\s*|\s$//g;
    my ($tmp) = $vic  =~ /([\w-]+)/;
    return("$vic", 1) unless $tmp;
    return("$tmp", 1) unless $vic eq $tmp;
    if($tmp =~ /^\d+$/) {
	#return("$tmp", 2) unless $dbh->selectrow_array("select 1 from primers where prmid= '$tmp'");
	return("$tmp", 2) unless $tmp <= @primr;
    } else {
	#$tmp = $dbh->selectrow_array("select prmid from primers where name = '$tmp'");
	return("$tmp", 2) unless $primr{$tmp};
    }
    return($tmp, 0);
}


sub chk_srcid {
    my ($vic) = $_[0]; $vic =~ s/^\s*|\s$//g;
    my ($tmp) = $vic  =~ /([\w\.:-]+)/;
    return($vic, 1) unless $tmp;
    return($tmp, 1) unless $vic eq $tmp;
    if($tmp =~ /^\d+$/) {
	return($tmp, 1) unless $tmp <= @srcid;
    } elsif($tmp =~ /^(\w+)[\s\.:-]+(\w+)$/) {
	return("$tmp", 3) unless $srcid{"$1-$2"};
    } elsif($tmp =~ /^[A-Za-z]+$/) {
	return("$tmp", 3) unless $srcid{"SE-$tmp"};
    }
    return($tmp, 0);
}


sub chk_sdate {
    my ($vic) = $_[0]; $vic =~ s/^\s*|\s$//g;
    my ($tmp) = $vic  =~ /(\d{4}-\d{2}-\d{2})/;
    return($vic, 1) unless $tmp && $vic eq $tmp;
    eval { Time::Piece->strptime($tmp, "%Y-%m-%d"); };
    return($tmp, 1) if $@;
    # err=3 if too old/future?
    return($tmp, 0);
}


sub feedback {
    my @err = @_; my $not = "";
    my $warn = '-- suspicious but OK';
    my $fail = '-- recheck or send ' . 
	a({href=>"mailto:$WADM?subject=$SITE: reg format"},'bug report');

    if ($err[0] == 1) { 
	$not = "Unreadable Run ID $fail";
    } elsif ($err[0] == 2) {
	# $not = "Duplicate Run ID $fail";
	return('Duplicate Run ID -- must be unique', (2,2,2,2,2));
    } elsif ($err[0] == 3) {
	$not = "Run ID not in standard format $warn";
    }

    if ($err[1] == 1) {
	$not .= br if $not; $not .= "Unreadable Sample ID $fail";
    } elsif ($err[1] == 3) {
	$not .= br if $not; $not .= "Sample ID not in standard format $warn";
    } elsif ($err[1] == 2) {
	$not .= br if $not; $not .= "Duplicate Sample ID $warn";
	$err[1] = 3;
    }


    if ($err[2] == 1) {
	$not .= br if $not; $not .= "Unreadable Primer ID $fail";
    } elsif ($err[2] == 2) {
	$not .= br if $not; $not .= "Unregistered Primer ID $fail";
	$err[2] = 1;
    }

    if ($err[3] == 1) {
	$not .= br if $not; $not .= "Unreadable Source ID $fail";
    } elsif ($err[3] == 3) {
	$not .= br if $not; $not .= "Unregistered Source ID $warn";
	#$err[3] = 1;
    }

    if ($err[4] == 1) {
	$not .= br if $not; $not .= "Unreadable Sample Date $fail";
    }

    ($not, @err);
}


sub Fret { 
    my $v = $_[0];
    return unless $v;
    print p(font({color=>'red'},$v));
}
