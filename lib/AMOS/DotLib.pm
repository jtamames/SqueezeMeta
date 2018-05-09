# $Id$
#
# DotLib.pm - set of procedures for generating .dot files
#

#  Copyright @ 2002, 2003, The Institute for Genomic Research (TIGR).  All
#  rights reserved.


=head1 Name

DotLib - library of routines for generating .dot files

=head1 Synopsis

    use DotLib;

=head1 Description

    A set of procedures used to create various .dot objects such as
file headers, file tails, components, nodes, edges, etc.

=cut

package AMOS::DotLib;

use strict;


BEGIN {
    use Exporter ();
    use vars qw(@EXPORT @EXPORT_OK @ISA %EXPORT_TAGS);

    @ISA         = qw(Exporter);
    @EXPORT      = qw(&printHeader
                      &printFooter
		      &printNode
		      &printEdge
		      &startCluster
		      &endCluster
		      );
    %EXPORT_TAGS = ();
    @EXPORT_OK   = ();
}

our $VERSION = '1.0'; 
our $REVISION = '$Revision$ ';
our $VERSION_STRING = "$VERSION ($REVISION)";

use vars @EXPORT;
use vars @EXPORT_OK;

=over 4

=item B<my $ret = printHeader($file, $type);>

Prints a .dot header for the type of output specified in the $type variable.
Allowable types are "printer", "plotter".  If $type is undefined or not
passed, it generates a default header.  Returns 1 upon successful 
completion and 'undef' otherwise.

Example:

    my $err = printHeader(\*STDOUT, "plotter");

=cut

sub printHeader
{
    my $file = shift;
    my $type = shift;

    print $file "digraph ROOT {\n";
    print $file "  rankdir = LR\n";
    print $file "  orientation = landscape\n";
    print $file "  ranksep = 0.3\n";
    print $file "  nodesep = 0.3\n";
    print $file "  fontsize = 8\n";
    print $file "  margin = \".2,.2\"\n";
	
    if ($type eq "printer"){
	print $file "  ratio = auto\n";
	print $file "  page = \"8.5,11\"\n";
    } elsif ($type eq "plotter"){
	print $file "  ratio = auto\n";
	print $file "  page = \"36,48\"\n";
    }
    
    print $file "\n";

    return 1;
} # printHeader


=item B<my $ret = printFooter($file);>

Prints a .dot footer (currently just a closed brace).  Returns 1 upon
successful completion and 'undef' otherwise.

Example:

    my $err = printFooter(\*STDOUT);

=cut

sub printFooter
{
    my $file = shift;

    print $file "}\n";
    
    return 1;
} # printFooter


=item B<my $ret = printNode($file, $id, $label, $ori);>

Prints a "contig" node with the specified id, label, and orientation.
If orientation is 1 then the node is a forward facing arrow, otherwise
it is a backward facing arror. Returns 1 upon successful completion
and 'undef' otherwise.

Example:

    my $err = printNode(\*STDOUT, $node_id, "$node_id ($node_len)", 1);

=cut   

sub printNode
{
    my $file = shift;
    my $id = shift;
    my $label = shift;
    my $ori = shift;
    my $angle;
    
    $id =~ s/(\W)/_/g;

    if ($ori == 1){
	$angle = -90;
    } else {
	$angle = 90;
    }

    print $file "    $id [ label = \"$label\" height = 0.2, fontsize = 8, shape = \"house\", orientation = $angle ]\n";

    return 1;

} # printNode


=item B<my $ret = printEdge($file, $nodeA, $nodeB, $label, $style);>

Prints an edge between two nodes with the specified label.  The style can
be any of the GraphViz acceptable styles ("dotted", "solid", "dashed", 
"invis") or undefined in which case the default is used. Returns 1 upon
successful completion and 'undef' otherwise.

Example:

    my $err = printEdge(\*STDOUT, $nodeA, $nodeB, "A to B", "invis");

=cut   

sub printEdge
{
    my $file = shift;
    my $nodeA = shift;
    my $nodeB = shift;
    my $label = shift;
    my $instyle = shift;
    my $style;

    $nodeA =~ s/(\W)/_/g;
    $nodeB =~ s/(\W)/_/g;

    if (defined $instyle){
	$style = "style = \"" . $instyle . "\"";
	if ($instyle eq "invis"){
	    $style .= " color = \"white\" ";
	}
    }

    print $file "    $nodeA -> $nodeB [ label =\"$label\" fontsize = 8 $style ]\n";

    return 1;
} # printEdge

=item B<my $err = startCluster($file, $id, $label);>

Starts a cluster in the .dot output file with the given label and id.
Returns 1 upon successful completion and 'undef' otherwise.

Example:
    
    my $err = startCluster(\*STDOUT, $clust_id, "first cluster");

=cut

sub startCluster
{
    my $file = shift;
    my $id = shift;
    my $label = shift;

    $id =~ s/(\W)/_/g;

    print $file "  subgraph cluster_$id {\n";
    print $file "    label = \"$label\"\n";

    return 1;
} # startCluster

=item B<my $err = endCluster($file);>

Ends a cluster in the .dot output.  Returns 1 upon successful
completion and 'undef' otherwise.

Example: 

    my $err = endCluster(\*STDOUT);

=cut

sub endCluster
{
    my $file = shift;

    print $file "  }\n";
   
    return 1;
} # endCluster


1;



