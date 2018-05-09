# $Id$
#

=head1 NAME

AMOS::Foundation - A library of perl routines for handling command line 
parameters and logging in a uniform fashion.

=head1 SYNOPSIS
    
use AMOS::Foundation;

$base = new AMOS::Foundation;

=head1 DESCRIPTION

AMOS::Foundation provides a set of routines useful in handling command
line parameters and logging.  These routines are meant to provide a uniform
house-keeping mechanism for AMOS related Perl scripts.

=cut

package AMOS::AmosFoundation;

use strict;
use Getopt::Long;
use POSIX;
use IO::File;

## configuration management            
our $VERSION = "1.00";
our $VERSION_STRING = "$VERSION (Build " . (qw/$Revision$/ )[1] . ")";
our @DEPEND = ();

BEGIN {
    use Exporter ();
    use vars qw(@EXPORT @EXPORT_OK @ISA %EXPORT_TAGS);

    require 5.006_00;                       # error if using Perl < v5.6.0  

    @ISA         = qw(Exporter);
    @EXPORT      = qw(&getOptions
                      &setHelpText
		      &setUsage
		      &setVersion
		      &getVersion
                      &setLogFile
                      &setLogLevel
                      &getHelpText
		      &getUsage
		      &getLogFile
		      &getLogLevel
		      &log
		      &bail
		      &new
                      &setMihaiMode
                      $VERSION
                      $VERSION_STRING
                      @DEPEND);

    %EXPORT_TAGS = ();
    @EXPORT_OK   = ();
}

use vars @EXPORT;
use vars @EXPORT_OK;


##############################################

=over

=item B<$base = new AMOS::Foundation; >

Creates a new AMOS::Foundation object.

=cut
sub new
{
    my $self = {};
    my $pkg = shift;

    # create the object
    bless $self, $pkg;
    
    $self->{"helpText"} = "";
    $self->{"Usage"} = "";
    $self->{"logLevel"} = 0;
    $self->{"logFile"} = "";
    $self->{"logStream"} = new IO::File;
    $self->{"logStream"}->open(">&STDERR") || die ("Cannot open STDERR: $!\n");

    return $self;
}


=item B<$base->setHelpText("this is some help info");>

Assigns the default text output when program is run with option -h or -help.

=cut

sub setHelpText
{
    my $self = shift;
    if (@_) {
	$self->{'helpText'} = shift;
    }

    return $self->{'helpText'};
}# setHelpText

=item B<$help = $base->getHelpText(); >

Returns the default text output when program is run with option -h or -help.

=cut

sub getHelpText
{
    my $self = shift;

    return $self->{'helpText'};
}# getHelpText

=item B<$base->setUsage("this is some help info"); >

Assigns the default text output when program is run with incorrect parameters.

=cut

sub setUsage
{
    my $self = shift;
    if (@_) {
	$self->{'Usage'} = shift;
    }

    return $self->{'Usage'};
}# setUsage

=item B<$usage = $base->getUsageText(); >

Returns the default text output when program is run with incorrect parameters.

=cut

sub getUsage
{
    my $self = shift;

    return $self->{'Usage'};
}# getUsage



=item B<$base->setLogFile("file.log"); >

Assigns the default location for the log file.  Log messages will be appended
to this file.

=cut

sub setLogFile
{
    my $self = shift;

    if ($self->{'logFile'} ne "" && exists $self->{'logStream'}){
	$self->{'logStream'}->close;  # close any previously opened file
    }
    
    if (@_) {
	$self->{'logFile'} = shift;
    }

    $self->{'logStream'}->open(">$self->{'logFile'}") || 
	die ("Cannot open $self->{'logFile'}: $!\n");

    return $self->{'logFile'};
}# setLogFile

=item B<$logfile = $base->getLogFile(); >

Retrieves location of the log file.

=cut

sub getLogFile
{
    my $self = shift;

    return $self->{'logFile'};
}# getLogfile

=item B<$base->setVersion("version"); >

Assigns the version to be returned when program called with -V option.

=cut

sub setVersion
{
    my $self = shift;

    if (@_) {
	$self->{'Version'} = shift;
    }

    return $self->{'Version'};
}# setVersion

=item B<$version = $base->getVersion(); >

Retrieves program version.

=cut

sub getVersion
{
    my $self = shift;

    return $self->{'Version'};
}# getVersion

=item B<$base->setLogLevel(2); >

Assigns the default log level.  Only messages at this level or lower will
be output.

=cut

sub setLogLevel
{
    my $self = shift;
    if (@_) {
	$self->{'logLevel'} = shift;
    }

    return $self->{'logLevel'};
}# setLogLevel

=item B<$logLevel = $base->getLogLevel(); >

Retrieves the current log level.

=cut

sub getLogLevel
{
    my $self = shift;

    return $self->{'logLevel'};
}# getLogLevel

=item B<$base->bail("message");>

Gracefully exit

=cut

sub bail
{
    my $self = shift;
    my $msg = shift;  # message to be sent to screen

    $self->log("ERROR: $msg", 0);

    if ($self->{'logFile'} ne ""){
	$self->{'logStream'}->close;
    }

    die();
}# die

=item B<$base->log("message", 5); >

Log a message to the log file if the default log level is above 5.

=cut

sub log
{
    my $self = shift;
    my $message = shift;
    my $level = shift;

    $level = 1 if (! defined $level);

    if ($self->{'logLevel'} < $level) {
	return; # only log at the appropriate level
    }

    my $date = strftime "%F %T| ", localtime;
    
    $self->{'logStream'}->print($date . $message . "\n");

    if ( defined $self->{'MihaiMode'} ) {
        print STDERR ($date . $message . "\n");
    }
}# die


=item B<$err = $base->getOptions("i=s" => \$in); >

Interface to Getopt::Long that provides special processing to -h, -V, -verbose.

=cut

sub getOptions
{
    my $self = shift;
    my @user_options = @_;
    my $logfile;
    my $verbose;
    my $version;
    my $help;
    my @AMOS_options = (
         "logfile=s" => \$logfile,
         "verbose=i" => \$verbose,
         "version|V" => \$version,
         "help|h" => \$help,
      );

    Getopt::Long::Configure('no_ignore_case');
    my $getopt_code = eval 'GetOptions (@user_options, @AMOS_options)';
    
    # In case of failure output usage
    if (! $getopt_code){
	print STDERR $self->{'Usage'}, "\n";
	exit(1);
    }

    if (defined $help){
	print $self->{'helpText'}, "\n";
	exit(0);
    }

    if (defined $version){
	print $self->{'Version'}, "\n";
	exit(0);
    }

    if (defined $verbose){
	$self->{'logLevel'} = $verbose;
    }

    if (defined $logfile){
	$self->{'logStream'}->close;
	$self->{'logFile'} = $logfile;
	$self->{'logStream'}->open(">$logfile") || 
	    die ("Cannot open $logfile: $!\n");
    }
    
    return $getopt_code;
}

=item B<$base->setMihaiMode(); >

Display log messages on stderr.

=cut

sub setMihaiMode
{
    my $self = shift;
    $self->{'MihaiMode'} = 1;
}

=back

=cut

END{}

1;
