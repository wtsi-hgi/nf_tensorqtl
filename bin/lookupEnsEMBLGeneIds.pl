#!/usr/bin/env perl

# look up gene annotations for a set of EnsEMBL IDs.
# using EnsEMBL mirror at the Sanger

##BEGIN {
##    $main::ensembl_version=81; # current ensembl version
##    unshift(@INC, "/nfs/users/nfs_h/hp3/sw/src/ensembl$main::ensembl_version/modules");
##}

use strict;
use warnings;

use Bio::EnsEMBL::Registry;

$main::Species ='homo sapiens';
# $main::ASSEMBLY_VERSION ='GRCh38';
# $main::ASSEMBLY_VERSION ='GRCh37';
$main::NameDB = "HGNC";

# $main::Species ='mus musculus';
# $main::ASSEMBLY_VERSION ='GRCm38';
# $main::NameDB = "MGI";

# command line arguments
my $genelist_file = shift;
my $oufnam = shift;
my $logfnam = shift;

(defined $genelist_file && defined $oufnam && defined $logfnam) ||
    die "call with <list of EnsEMBL IDs> <ouput file name> <outfil failed Ids>\n";

# load ENSEMBL registry
my $registry = "Bio::EnsEMBL::Registry";

$registry->load_registry_from_db(
  -host => 'ensembldb-mirror.internal.sanger.ac.uk',
  -user => 'anonymous'
);

my $gene_adaptor = $registry->get_adaptor( $main::Species, 'Core', 'Gene' );
$gene_adaptor or die "Can't get gene adaptors\n";

##$main::slice_adaptor = $registry->get_adaptor( $main::Species, 'Core', 'Slice' );
##$main::slice_adaptor or die "Can't get slice adaptors\n";

# load list of gene names
open (my $infh, "<", $genelist_file) or die "could not open file $genelist_file\n";
open (my $oufh, ">", $oufnam) or die "could not open file $oufnam\n";
open (my $logfh, ">", $logfnam) or die "could not open file $logfnam\n";
print $oufh "query_id\tensembl_id\tdisplay_id\tchromosome\tstart\tend\tstrand\thgnc_symbols\n";
GENE: while (<$infh>) {
    my @fld = split;
    my $genenam = uc $fld[0];
    if (substr($genenam,0,1) eq '#') {
	     next GENE;
    }
    ##print "lookup gene: $genenam\n";
    my $gene = $gene_adaptor->fetch_by_stable_id($genenam);
    if (!defined $gene) {
	     $gene = $gene_adaptor->fetch_by_display_label($genenam);
    }
    if (!defined $gene) {
	     print "gene $genenam not found\n";
       print $logfh "$genenam\n";
	     next GENE;
    }
    #my $coosysnam = $gene->coord_system_name();
    my $ensid = $gene->stable_id();
    my $gen_id = $gene->display_id();
    #print "Assembly version: $coosysnam\n";
    my $chrnam = $gene->seq_region_name();
    my $start = $gene->start();
    my $end = $gene->end();
    my $strand = $gene->strand();
    ##my $xrefs = $gene->get_all_DBEntries("HGNC");
    my $xrefs = $gene->get_all_DBEntries();
    my ($gen_symbol, $synonyms) = getHGNCfromDBEntries($xrefs);
    #print $oufh "$ensid\t$gen_symbol";
    if ($chr =~ /^[1-9,X,Y]$/ || $chr =~ /^1[1-9]/ || $chr =~ /^2[12]/ || $chr == 'MT') {
        $chr = 'chr'.$chr
    }
    print $oufh "$genenam\t$ensid\t$gen_id\t$chrnam\t$start\t$end\t$strand\t$gen_symbol";
    foreach my $id (@{$synonyms}) {
	     print $oufh "\t$id";
    }
    print $oufh "\n";
}

close $logfh;
close $oufh;
close $infh;

exit 0;

sub getHGNCfromDBEntries
{
    my $db_entries = shift;
    my $symbol = "N/A";
    my $synonyms;
ENTRY:
    foreach my $dbe ( @{$db_entries} ) {
	#printf "\tXREF %s (%s)\n", $dbe->display_id(), $dbe->dbname();
	if ($dbe->dbname() eq $main::NameDB) {
	    $symbol = $dbe->display_id();
	    $synonyms = $dbe->get_all_synonyms();
	    last ENTRY;
	}
    }
    return ($symbol, $synonyms)
}
