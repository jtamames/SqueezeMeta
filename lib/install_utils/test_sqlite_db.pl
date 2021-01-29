use DBI;

$databasepath = $ARGV[0];
$lca_db    = "$databasepath/LCA_tax/taxid.db";

my $dbh = DBI->connect("dbi:SQLite:dbname=$lca_db","","",{ RaiseError => 1, sqlite_open_flags => SQLITE_OPEN_READONLY }) or die $DBI::errstr;
$dbh->sqlite_busy_timeout( 120 * 1000 );
my $query="select * from taxid LIMIT 1";
my $sth = $dbh->prepare($query);  
$sth->execute();
my @list;
while(@list=$sth->fetchrow()) { } #{ print "@list\n"; }
