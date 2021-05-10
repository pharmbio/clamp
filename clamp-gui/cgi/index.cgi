#!/usr/bin/perl -T
#
# Name: Philadelphia portal
# Auth: Wesley Schaal, wesley.schaal@farmbio.uu.se
# When: 2014-10-09
use lib ('../perl/lib/perl5','.');
use DBI; use strict; use warnings; use v5.10;
use CGI qw(:all *table *Tr *td *blockquote -nosticky);
use CGI::Carp qw(fatalsToBrowser);  autoEscape(undef);
use List::Util qw(first);
use Time::Piece; use Time::Seconds;
use common; # local hacks

my $appl = lc $SITE; # common
my $desc = $DESC;    # common
my $data = 'uploads';
my $D=0;

my $dbf = "../data/${appl}.sqlite";
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbf",'','');

# my $db   = $appl; 
# # my $grp  = server_name() eq 'pele.farmbio.uu.se' ? "${appl}loc" : "${appl}bmc";
# my $grp  = "${appl}loc";
# my $dbh  = DBI->connect("DBI:mysql:$db"
# 	   . ";mysql_read_default_file=../data/.my.cnf"
# 	   . ";mysql_read_default_group=$grp") ||
#            die "Can't connect to $db: $DBI::errstr";

my ($min,$hour,$day,$mon,$year) = (localtime)[1..5];
$year += 1900; $mon++;
my $today = sprintf("%4d-%02d-%02d", $year,$mon,$day);
my $now   = sprintf("Fetched: %4d-%02d-%02d at %02d:%02d",
		 $year,$mon,$day,$hour,$min);

#for my $v (qw(rec pat submit))  { $$v = param($v); }
my $rec = param('rec'); my $pat = param('pat'); my $all = param('all');
my $submit = param('submit');
ShowHead();
if   ($rec)      { ShowReg(); }
#elsif($pat)      { ShowPat(); }
elsif($all)      { ShowAll(); }
else             { FindHit(); }
ShowFoot();


sub FindHit {
    print start_blockquote;
    print h2("Overview of $desc");
    print p(font({color=>'#D05050'},'Enter search terms into the the form to ' .
		 'retreive data. Empty fields means "anything" for that field. '.
		 '')) unless $submit;

    print start_form;
    print start_table({border=>1,align=>'center',width=>'95%',
		       cellspacing=>0,cellpadding=>'10'});
    print start_Tr, start_td;
    print start_table({width=>'95%',cellpadding=>'5'}), start_Tr;
    print td(['Sample ID' . br .
	      textfield(-name=>'samid',size=>'12',maxlength=>'12'),
	      'Run ID' . br .
	      textfield(-name=>'runid',size=>'16',maxlength=>'64'),
	      'Date (yyyy-mm-dd)' . br .
	      textfield(-name=>'sdate',size=>'18',maxlength=>'32'),
	      ]);
    print end_Tr, end_table, end_td, end_Tr, start_Tr, start_td;
    print start_table({width=>'95%',cellpadding=>'5'}), start_Tr;
    print td(['Mutation' . br .
	      textfield(-name=>'mutat',size=>'18',maxlength=>'32'),
	      'Reference' . br .
	      textfield(-name=>'ref',size=>'18',maxlength=>'32'),
	      'Uploaded' . br .
	      popup_menu(-name=>'seq',-values=>['-','Yes','No'],-default=>'Yes'),
	      ]);
    print end_Tr, end_table;
    print end_td, end_Tr, end_table;
    print table({align=>'center',width=>'50%',cellspacing=>'10'},
		Tr(th([submit(-name=>'submit',value=>'Search'),
		       defaults('Reset')])));
    print end_form, br;

    PrepSQL() if $submit;

    print end_blockquote;
}


sub PrepSQL {
    my @vars = qw(samid runid sdate seq mutat ref);
    print hr({width=>'65%'});

    my $where = '';
    for my $v (@vars)  {
	my $p = param($v);
	$p =~ s/^\s*|\s*$//g; $p =~ s/^-$//;
	if ($v eq 'sdate' && $p) {
	    # todo: accept 'yyyy' -> 'yyyy-01-01' (with modifiers)
	    $p =~ y/0-9<>-//cd; next unless $p;
	    $p =~ s/^-/</; $p =~ s/(.+)-$/>$1/; 
	    if     ($p =~ /^(\d{4})-?(\d\d)-?(\d\d)$/) { 
		$where .= "$v  = '$1-$2-$3' and ";
	    } elsif($p =~ /^(\d{2})-?(\d\d)-?(\d\d)$/) { 
		$where .= "$v  = '20$1-$2-$3' and ";
	    } elsif($p =~ /^([<>])\s*(\d{4})-?(\d\d)-?(\d\d)$/) { 
		$where .= "$v $1 '$2-$3-$4' and ";
	    } elsif($p =~ /^\s*([<>])\s*(\d{2})-?(\d\d)-?(\d\d)$/) { 
		$where .= "$v $1 '20$2-$3-$4' and ";
	    } elsif($p =~ /^(\d{4})-?(\d\d)-?(\d\d)[\s-]+(\d{4})-?(\d\d)-?(\d\d)$/) {
		$where .= "($v between '$1-$2-$3' and '$4-$5-$6'";
	    } elsif($p =~ /^(\d{2})-?(\d\d)-?(\d\d)[\s-]+(\d{2})-?(\d\d)-?(\d\d)$/) {
		$where .= "($v between '20$1-$2-$3' and '20$4-$5-$6') and ";
	    } else { Fret("$v = $p"); }
	} elsif ($v eq 'seq' && $p) {
	    my $subq = '(select distinct regid from files where ftype = "txt")';
	    if    ($p eq 'Yes') { $where .= "regid in $subq and "; }
	    elsif ($p eq 'No')  { $where .= "regid not in $subq and "; }
	    else { Fret($v); }
	} elsif ($v eq 'mutat' && $p) {
	    my($q, $sub) = ('mutation','');
	    $p =~ y/A-Za-z0-9\!\*-//cd; next unless $p; $p = uc $p;
	    $p =~ s/unk\w*/rs/i; $p =~ s/rs/rs*/i;
	    if      ($p =~ s/^[\!-](\w*\*)/$1/) { 
		$p =~ y/*/%/;
		$sub .= "$q like '$p' and detection = 'negative'";
	    } elsif ($p =~ /\*/) { 
		$p =~ y/*/%/;
		$sub .= "$q like '$p' and detection = 'positive'";
	    } elsif ($p =~ /^[\!-](\w+)/) { 
		$sub .= "$q = '$1' and detection = 'negative'";
	    } elsif ($p =~ /^\w+/) { 
		$sub .= "$q  = '$p' and detection = 'positive'";
	    } else { Fret("$v = $p"); }
	    if($sub) {
		#$sub =~ s/and\s*$//;
		$where .= "regid in (select regid from results where $sub) and ";
	    }
	} elsif ($v eq 'ref' && $p) {
	    $p =~ y/0-9a-zA-Z_\*\!-//cd;
	    $p =~ y/*/%/;
	    $where .= "regid in (select regid from logs where $v like '$p') and ";
	} elsif ($p) {
	    $p =~ y/0-9a-zA-Z_\*\!-//cd;
	    $p =~ y/*/%/;
	    $where .= "$v like '$p' and ";
	# } elsif ($v ~~ [qw(samid runid)] && $p) {
	#     $p =~ y/0-9<>,!-//cd;
	#     # $p =~ y/0-9a-zA-Z_<>,!-//cd;
	#     next unless $p;
	#     my $q = $v eq 'foobar' ? "a.$v" : $v;
	#     if   ($p =~ /^\d+$/)               { $where .= "$q  = $p and "; }
	#     elsif($p =~ /^[\!-]\s*(\d+)$/)     { $where .= "$q != $1 and "; }
	#     elsif($p =~ /^\<\s*(\d+)$/)        { $where .= "$q  < $1 and "; }
	#     elsif($p =~ /^\>\s*(\d+)$/)        { $where .= "$q  > $1 and "; }
	#     elsif($p =~ /^(\d+)\s*-\s*(\d+)$/) { $where .= "($q between $1 and $2) and "; }
	#     elsif($p =~ /\d+,[\d\s,]*\d+/)     { $where .= "$q in ($p) and "; } # ",,"
	#     else { Fret("$v = $p"); }
	}
	print p("$v = '\%$p\%'") if $D>1;
    }

    if($where) {
	$where = "where $where";
	$where =~ s/and\s*$//;
    }

    my $sql = "select distinct regid from samples $where"; print p($sql) if $D;
    my $hits = $dbh->selectcol_arrayref($sql); 
    my @hits = defined $hits ? @$hits : ();
    print p(@hits) if $D;     #return(1);

    #@hits ? ListRec(@hits) : print p("No matches found");
    @hits ? ListAll(@hits) : print p("No matches found");
}


# sub ListRec {
#     my $recs = join ',', grep { /^\d+$/ } @_;
#     my $sql = "select regid,samid,runid,sdate,(select group_concat(mut) from ".
# 	"mutants where eid=regid) mutat from samples where regid in ($recs)";
#     my $sth = $dbh->prepare($sql); $sth->execute(); print(p($sql)) if $D;

#     print start_blockquote, br, start_form;
#     print submit(-name=>'all',value=>'Details');
#     print start_table({cellspacing=>0,cellpadding=>5,border=>1});
#     print Tr(th([qw(DB Patient Sample Date Mutations Other)]));
#     while(my ($regid,$samid,$runid,$sdate,$mutat) = $sth->fetchrow_array) {
# 	$mutat =~ y/,/ /;
# 	print Tr(td({align=>'center'}, submit(-name=>'rec',value=>$regid)) .
# 		 td({align=>'center'}, submit(-name=>'pat',value=>"D:$samid")) .
# 		 td("R$runid") . td($sdate) . td($mutat) . 
# 		 td({align=>'center'},chr(int(rand(26))+65)));
#     }
#     print end_table, p(font({size=>'-1',color=>'gray'},$now));
#     print hidden('recs',$recs);
#     print end_form, end_blockquote;
# }


sub ListAll {
    my $recs = join ',', grep { /^\d+$/ } @_; print p("Recs: $recs") if $D;
    # my $recs = param('recs'); print "HITS: $recs" if $D;
    Fret('Internal Error 21') unless $recs =~ /^\d+(?:,\d+)*$/;

    # my $sql = "select regid,mutation as mut,round(100*freq,1) as pct " .
    #  	"from results where regid in ($recs) and detection = 'positive'";
    # my $sql = "select regid, case when instr(mutation, 'rs') then " .
    # 	"'Unknown' else mutation end as mut, round(100*freq,1) as pct " .
    # 	"from results where regid in ($recs) and detection = 'positive'";
    my $sql = "select regid,mutation as mut,round(100*freq,1) as pct " .
	"from results where regid in ($recs) and detection = 'positive' " .
	'and mutation not like "rs%"';
    print p($sql) if $D; # use 'case' for compatability (sqlite,mysql)
    my $res = $dbh->selectall_hashref($sql, [ qw(regid mut) ]);

    my @mut = uniq(map { keys %{$res->{$_}} } %$res);
    @mut = map {$_->[0]} sort {$a->[1]<=>$b->[1]} map {[$_, $_=~/(\d+)/]} @mut;
    print p(join(',',@mut)) if $D;

    $sql = "select regid, count(*) as cnt from results where " .
	"detection = 'unresolved' and regid in ($recs) group by regid";
    my $unr = $dbh->selectall_hashref($sql, 'regid');

    $sql = 'select regid, count(*) as cnt from results where ' .
	'mutation like "rs%" and regid in ' . "($recs) group by regid";
    my $unk = $dbh->selectall_hashref($sql, 'regid');

    $sql = "select samid,regid,runid,sdate,(select 1 from files b where " .
    	"a.regid=b.regid and ftype = 'txt' order by stamp desc limit 1) " .
    	"uploaded from samples a where regid in ($recs)"; print p($sql) if $D;
    my $sth = $dbh->prepare($sql); $sth->execute();

    $sql = "select regid, ref from logs";
    my $ref = $dbh->selectall_hashref($sql, 'regid');

    print start_blockquote, br, start_form;
    print hidden('recs', $recs);
    print start_table({cellspacing=>0,cellpadding=>5,border=>1});
    print Tr(th({valign=>'top'},['Details', 'Sample ID', 'Run ID', 'Reference',
		 'Unresolved', @mut, 'Date']));
    while(my ($samid,$regid,$runid,$sdate,$uploaded) = $sth->fetchrow_array) {
	print start_Tr, td({align=>'center'},
			   [submit(-name=>'rec',value=>$regid),$samid,$runid]);

	print td({align=>'center'},$ref->{$regid}{ref});
	
	# my $unr_row = $unr->{$regid}{cnt}; my $unr_clr;
	# my $unr_row = $unr->{$regid}{cnt} ? "QC failed" : '&nbsp;'; 
	my $unr_row = $unr->{$regid}{cnt} ? "QC failed" : ''; 
	my $unr_clr;
	if(!$uploaded)  { $unr_clr = '#d0d0d0'; }
	elsif($unr_row) { $unr_clr = '#ffff99'; }
	else            { $unr_clr = '#ffffff'; }
	print td({align=>'center',bgcolor=>$unr_clr},$unr_row);

	# my $unk_row = $unk->{$regid}{cnt}; my $unk_clr;
	# if(!$uploaded)  { $unk_clr = '#d0d0d0'; }
	# elsif($unk_row) { $unk_clr = '#99ff99'; }
	# else            { $unk_clr = '#ffffff'; }
	# print td({align=>'center',bgcolor=>$unk_clr},$unk_row);

	for my $mut (@mut) {
	    my $val = $res->{$regid}{$mut}{pct};
	    my $clr = sprintf "#ff%02xff", int(255-2.55*$val);
	    $clr = '#d0d0d0' unless $uploaded;
	    print td({align=>'center',bgcolor=>$clr},$val);
	}
	print td({align=>'center'},[$sdate]), end_Tr;
    }
    print end_table, end_form;
    # print p(font({color=>"#909090"}, url(-path_info=>1,-query=>1)));
    print p(font({color=>"#909090"}, url(-query=>1)));
    print end_blockquote;
}

sub ShowReg {
    $rec =~ y/0-9//cd; unless ($rec) { FindHit(); ShowFoot(); exit(0); }
    my $recs = param('recs');
    my $sql = "select mutation,sequence,wt_reads,mut_reads,other_reads,freq,".
	"detection from results where regid = $rec order by sortby"; print p($sql) if $D;
    my $sth = $dbh->prepare($sql); $sth->execute(); my @cols = @{ $sth->{NAME} };

    $sql = "select samid,runid,sdate from samples where regid = $rec";
    my ($samid,$runid,$sdate) = $dbh->selectrow_array($sql);

    # $sql = "select max(case when ftype = 'txt' then versn end) as txt, " .
    # 	"max(case when ftype = 'seq' then versn end) as seq, max(case when " .
    # 	"ftype = 'pdf' then versn end) as pdf from files where regid = $rec";
    # my ($txt,$seq,$pdf) = $dbh->selectrow_array($sql);

    # my $stf = $dbh->prepare('select orig || "." || SUBSTR("00" || versn, -2, 2) '.
    # 			    'from files where regid = ? and ftype = ? '.
    # 			    'group by regid having max(versn)');

    # my $stf = $dbh->prepare('select runid, orig, versn from files join samples '.
    # 			    'using (regid) where regid = ? and ftype = ? group '.
    # 			    'by regid having max(versn)');

    my $stf = $dbh->prepare('select runid, orig from files join samples using '.
			    '(regid) where regid = ? and ftype = ? group by '.
			    'regid having max(versn)');

    my @txt = $dbh->selectrow_array($stf,undef,($rec,'txt'));
    my @seq = $dbh->selectrow_array($stf,undef,($rec,'seq'));
    my @pdf = $dbh->selectrow_array($stf,undef,($rec,'pdf'));
    my @dst = $dbh->selectrow_array($stf,undef,($rec,'dist'));
    my @dtt = $dbh->selectrow_array($stf,undef,($rec,'dtab'));
    my @log = $dbh->selectrow_array($stf,undef,($rec,'log'));

    print start_blockquote, br;
    print p("Rec: $rec"), p("Recs: $recs") if $D;
    Nav($rec, $recs);

    print start_table({width=>'50%',align=>'center',bgcolor=>'#F0F0F6'});
    print start_Tr, start_td;
    print start_table({align=>'center',width=>'90%'});
    print Tr(th(["Sample ID","Run ID", "Date"]));
    print Tr(td({align=>'center'},[$samid,$runid,$sdate]));
    print end_table, end_td, end_Tr, end_table, br, br; #, hr({width=>'65%'});

    my $fn = sprintf "%s/%s-%06d", $data, $appl, $rec; print p($fn) if $D;
    print start_table({width=>'50%',align=>'center'});
    print Tr(th({colspan=>4},'Downloads:')), start_Tr;
    print $txt[1] ? 
#	td(a({href=>sprintf("%s-%02d.%s",$fn,$txt,'txt')},'Results')) : 
#	td(a({href=>sprintf("%s/%s/%s.%02d",$data,@txt)},'Results')) : 
	td(a({href=>sprintf("%s/%s/%s",$data,@txt)},'Results')) : 
	th(font({color=>'#C0C0C0'},'Results'));
    print $seq[1] ? 
#	td(a({href=>sprintf("%s-%02d.%s",$fn,$seq,'fastq.gz')},'Sequence')) : 
#	td(a({href=>sprintf("%s/%s/%s.%02d",$data,@seq)},'Sequence')) : 
	td(a({href=>sprintf("%s/%s/%s",$data,@seq)},'Sequence')) : 
	th(font({color=>'#C0C0C0'},'Sequence'));
    print $pdf[1] ? 
#	td(a({href=>sprintf("%s-%02d.%s",$fn,$pdf,'pdf')},'Details')) : 
#	td(a({href=>sprintf("%s/%s/%s.%02d",$data,@pdf)},'Details')) : 
	td(a({href=>sprintf("%s/%s/%s",$data,@pdf)},'Coverage')) : 
	th(font({color=>'#C0C0C0'},'Coverage'));
    print $dtt[1] ? 
	td(a({href=>sprintf("%s/%s/%s",$data,@dtt)},'Clonal txt')) : 
	th(font({color=>'#C0C0C0'},'Clonal txt'));
    print $dst[1] ? 
	td(a({href=>sprintf("%s/%s/%s",$data,@dst)},'Clonal pdf')) : 
	th(font({color=>'#C0C0C0'},'Clonal pdf'));
    print $log[1] ? 
	td(a({href=>sprintf("%s/%s/%s",$data,@log)},'Log')) : 
	th(font({color=>'#C0C0C0'},'Log'));    
    print end_Tr, end_table, hr({width=>'65%'}), br;

    print start_table({align=>'center',cellspacing=>0,cellpadding=>5,border=>1});
    print Tr(th(\@cols));
    while(my @row = $sth->fetchrow_array) {
	print Tr(td(\@row));
    }
    print end_table;
    print p(font({color=>"#909090"},'This sample: ' .br.url() . "?rec=$rec"));
    print end_blockquote;
}


sub Nav {
    use List::Util qw(first);
    my ($n, $a) = @_; my @v = split(/,/,$a);
    $n = first { $v[$_] == $n } 0..$#v;
    # ($n) = grep { $v[$_] ~~ $n } 0 .. $#v;

    my ($typ,$tag) = ('rec','Sample');
    my $prev = $n        ? $v[$n-1] : '';
    my $next = $n == $#v ? '' : $v[$n+1];

    print start_form, hr({width=>'65%'});
    print start_table({align=>'center'}),start_Tr;
    print start_td({bgcolor=>'#f0e0ff'});
    print $prev ? submit(-name=>$typ,value=>$prev) :
        font({color=>'#606060'},' First ');
    print ' &nbsp; &nbsp; ';
    print b(" $tag $v[$n] ");
    print ' &nbsp; &nbsp; ';
    print $next ? submit(-name=>$typ,value=>$next) :
        font({color=>'#606060'},' Last ');
    print end_td,td(' &nbsp; &nbsp; ');
    print td(' &nbsp; &nbsp; ');
    print td({bgcolor=>'#e0f0ff'},'&nbsp;' .
	     a({href=>url},'New Search'));
    # submit(-name=>'list',value=>'Search results')); # restore srch params
    print end_Tr,end_table,hr({width=>'65%'});
    print hidden('recs',join(',',@v));
    print end_form, br;
}


sub Fret {
    my $v = $_[0];
    return unless $v;
    $v = "Ignoring unexpected value '$v'";
    print font({size=>'+1',color=>'green'},$v);
}


sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}


# sub ShowPat {
#     #$pat =~ y/0-9//cd; 
#     unless ($pat) { FindHit(); ShowFoot(); exit(0); }
#     my $recs = param('recs'); #Fret('Internal Error 12') unless $recs =~ /^\d+(?:,\d+)*$/;
#     print p("pat= $pat".br."recs= $recs") if $D;

#     my @recs = grep { /^\d+$/ } split /,/, $recs;
#     my $sql = "select distinct samid from samples where regid in ($recs) order by samid";
#     print p($sql) if $D;
#     my @hits = @{ $dbh->selectcol_arrayref($sql) };

#     my $n = first {$hits[$_] eq "$pat"} 0..$#hits;
#     print "HITS: $n ; $#hits ; $pat ; ", join(',',@hits) if $D;

#     Nav('p', $n, @hits);

#     ####################
#      my $sql = "select count(*) n, sdate, ruser from samples " .
# 	 "where samid=$pat group by samid having min(sdate)";
#      my ($n,$sdate,$ruser) = $dbh->selectrow_array($sql);
    
#     print start_table({width=>'50%',align=>'center',bgcolor=>'#F0F0F6'});
#     print start_Tr, start_td;
#     print start_table({align=>'center',width=>'90%'});
#     print Tr(th([qw(Patient Samples Registered Registrar)]));
#     print Tr(td({align=>'center'},[sprintf("%05d",$pat),$n,$sdate,$ruser]));
#     print end_table, end_td, end_Tr, end_table, br, hr({width=>'65%'});
#     ####################

#     $sql = "select eid,mut,pct from mutants left join samples on (regid=eid) " .
# 	"where samid=$pat"; print p($sql) if $D;
#     my $res = $dbh->selectall_hashref($sql, [ qw(eid mut) ]); 
#     my @mut = uniq(map { keys %{$res->{$_}} } %$res);
#     @mut = map {$_->[0]} sort {$a->[1]<=>$b->[1]} map {[$_, $_=~/(\d+)/]} @mut;
#     print p(join(',',@mut)) if $D;

#     $sql = "select regid,runid,sdate,(select files from seqs b where " .
#     	"a.regid=b.regid order by stamp desc limit 1) seq " .
#     	"from samples a where samid=$pat"; print p($sql) if $D;
#     my $sth = $dbh->prepare($sql); $sth->execute();

#     print start_blockquote, br, start_form;
#     print start_table({cellspacing=>0,cellpadding=>5,border=>1});
#     print Tr(th(['DB', 'Sample', @mut, 'Seq', 'Date']));
#     while(my ($regid,$runid,$sdate,$seq) = $sth->fetchrow_array) {
#     	my @val = map { $res->{$regid}{$_}{pct} } @mut;
# 	print Tr(td({align=>'center'},[submit(-name=>'rec',value=>$regid),
# 				       "R$runid",@val,$seq,$sdate]));
#     }
#     print end_table, end_form, end_blockquote;
# }


