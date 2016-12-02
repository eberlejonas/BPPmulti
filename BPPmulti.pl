#!/usr/bin/env perl

use warnings;
use strict;

#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# A script for automatically running iBPP 2.1.2 and BPP3 multiple
# times
# 1. with pairwise combinations of X tau and Y theta priors
# 2. on multiple data sets (prior, molecular, trait, integrative)
# 3. on multiple guide tree topologies
# 4. in parallel
#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

use threads;

#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# USER SETTINGS
#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# set global variables
my $maxThreads = 4;                    # How many threads can run at one time

# executables
my $BPP    = './bpp';                  # obligatory
my $iBPP   = './ibpp';                 # obligatory

### working dir and input files
# subfolders with control files and (i)BPP-output will be created in the working directory
my $wdir          = '$HOME/BPPmulti/'; # use trailing slash
my $seqfile       = 'sequences.txt';   # all relative to working directory
my $Imapfile      = 'Imap.txt';
my $traitfile     = 'morph.txt';       # used for iBPP
my $heredityfile  = 'heredity.txt';
my $locusratefile = '';
my @trees         = ('(sp9,((sp11,spX2),((sp12,(sp1,sp2)),(sp6,sp10))));',
                     '(sp9,((sp11,spX2),(sp10,(sp6,(sp12,(sp1,sp2))))));',
                     '(sp9,((sp11,spX2),((sp10,(sp1,sp6)),(sp12,sp2))));'
                    );

#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# (i)BPP settings

my $speciesmodelprior       = 1; # 0: uniform labeled histories | 1: uniform rooted trees
my $speciesDelimitationBPP  = '1 1 1.5 1.0';
my $speciesDelimitationIBPP = '1 1 1.5 1.0 0 1';

my $cleandata = 0;  # remove sites with ambiguity data? 0|1 means off|on
my $nloci     = 4;  # number of data sets in seqfile
my $ntraits   = 4;  # number of trait variables
my $nindT     = 68; # total # individuals for which trait data is available

my @thetapriors = ('1 10', '1 20', '2 2000');       # for standard-prior testing: list of gamma(a, b) for theta
my @taupriors   = ('1 10 1', '1 20 1', '2 2000 1'); # for standard-prior testing: list of gamma(a, b) for root tau & Dirichlet(a) for other tau's
my $nu0         = 0;                                # nu and kappa: parameters for prior on traits
my $kappa0      = 0;                                # nu0 = 0 and kappa0 = 0 for non-informative prior

my $heredity   = '2';                               # ('0': No variation, '1': estimate, '2': from file) & a_gamma b_gamma (if 1)
my $locusrate  = '1 2.0';                           # ('0': No variation, '1': estimate, '2': from file) & a_Dirichlet (if 1)
my $sequenceerror;                                  # NOT IMPLEMENTED NOW sequencing errors: gamma(a, b) prior; no value for switching it of
my $finetune   = '1: 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01'; # auto on|off: GBtj, GBspr, theta, tau, mix, locusrate, seqerr, traitHsq

# MCMC settings: total number of generations = burnin + samplefreq * nsample
my $print    = 1;     # 0|1 means off|on
my $burnin   = 10000;
my $sampfreq = 2;
my $nsample  = 25000; # the authors stated in the manual (?) that 5e4 iterations should be sufficient for most real world data. Check ESS values in $wdir/final_output.txt!

# switch desired tasks on|off
my $UseBPP_prior    = 1; # use BPP without data (prior only)
my $UseBPP          = 1; # use BPP for molecular data alone
my $UseBPP_NNI      = 1; # use BPP for molecular data alone with species tree estimation
my $UseiBPP_intgr   = 1; # use iBPP for an integrative analysis
my $UseiBPP_trait   = 1; # use iBPP for an analysis on continuous traits only

my $executeAnalyses = 1; # set to 1 if you want to run the analyses
my $repeats         = 5; # set the number of repeats for each analysis

#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# USER SETTINGS ARE FINISHED HERE
#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~




# rudimentary user input checking

if (scalar @thetapriors != scalar @taupriors ) {
	print "WARNING: Number of alpha/beta value pairs in theta priors does not match those in tau priors!\n";
}


# check if necessary files exist
unless (-e "$wdir$seqfile") {
	die "Sequence-file '$wdir$seqfile' doesn't exist!\n"
}
unless (-e "$wdir$Imapfile") {
	die "Imap-file '$wdir$Imapfile' doesn't exist!\n"
}
unless (-e "$wdir$traitfile") {
	die "Trait-file '$wdir$traitfile' doesn't exist!\n"
}
unless (-e "$wdir$heredityfile") {
	die "Heredity-file '$wdir$heredityfile' doesn't exist!\n"
}
unless (-e "$wdir$locusratefile") {
	die "Locus-Rate-file '$wdir$locusratefile' doesn't exist!\n"
}


# starting BPP batch

print "

	" . scalar localtime() . ": Starting BPP batch processing
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~\n";
print "Changing to working directory\n";
chdir $wdir or die "Can't cd into working directiory '$wdir': $!\n";

# internal global variables
my @prior_means;
my @prior_parameters;
my @analyses;
my @threads;
my $speciesdelimitation;
my $speciestree;
my $useseqdata;
my $usetraitdata;


print "Parsing Imap file\n";
my $SpeciesNumberString = &parseImap ($Imapfile);


# write control files
for (my $rep_i = 0; $rep_i < $repeats; $rep_i++) {
	for (my $th_i = 0; $th_i < scalar @thetapriors; $th_i++) {
		for (my $tau_i = 0; $tau_i < scalar @taupriors; $tau_i++) {
			my $i = 1;
			foreach my $tree (@trees) {
				my $analysis_base;
				my $wdir_rep;
				$wdir_rep      = $wdir.'repeat'.$rep_i.'/';
				$analysis_base = sprintf "std.theta%d-tau%d.Tree%d", $th_i+1, $tau_i+1, $i;
				
				# prior only
				if ($UseBPP_prior == 1) {
					my $analysis = $analysis_base . '.BPP.prior';
					push (@analyses, $analysis) if $rep_i == 0;
					$speciesdelimitation = $speciesDelimitationBPP;
					$speciestree         = 0;
					$useseqdata          = 0;
					&makeControlBPP ($wdir_rep, $analysis, $SpeciesNumberString, $tree, $thetapriors[$th_i], $taupriors[$tau_i]);
				}
				
				# molecular data
				if ($UseBPP == 1) {
					my $analysis = $analysis_base . '.BPP.mol';
					push (@analyses, $analysis) if $rep_i == 0;
					$speciesdelimitation = $speciesDelimitationBPP;
					$speciestree         = 0;
					$useseqdata          = 1;
					&makeControlBPP ($wdir_rep, $analysis, $SpeciesNumberString, $tree, $thetapriors[$th_i], $taupriors[$tau_i]);
				}
				
				# molecular data + species tree
				if ($UseBPP_NNI == 1) {
					my $analysis = $analysis_base . '.BPP.mol-NNI';
					push (@analyses, $analysis) if $rep_i == 0;
					$speciesdelimitation = $speciesDelimitationBPP;
					$speciestree         = 1;
					$useseqdata          = 1;
					&makeControlBPP ($wdir_rep, $analysis, $SpeciesNumberString, $tree, $thetapriors[$th_i], $taupriors[$tau_i]);
				}
				
				# trait data
				if ($UseiBPP_trait == 1) {
					my $analysis = $analysis_base . '.iBPP.trait';
					push (@analyses, $analysis) if $rep_i == 0;
					$speciesdelimitation = $speciesDelimitationIBPP;
					$useseqdata          = 0;
					$usetraitdata        = 1;
					&makeControliBPP ($wdir_rep, $analysis, $SpeciesNumberString, $tree, $thetapriors[$th_i], $taupriors[$tau_i]);
				}
				
				# integrative
				if ($UseiBPP_intgr == 1) {
					my $analysis = $analysis_base . '.iBPP.intgr';
					push (@analyses, $analysis) if $rep_i == 0;
					$speciesdelimitation = $speciesDelimitationIBPP;
					$useseqdata          = 1;
					$usetraitdata        = 1;
					&makeControliBPP ($wdir_rep, $analysis, $SpeciesNumberString, $tree, $thetapriors[$th_i], $taupriors[$tau_i]);
				}
				$i++;
			}
		}
	}
}

# execute standard prior runs
if ($executeAnalyses == 1) {
	print "

	" . scalar localtime() . ": Starting execution of standard prior runs
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~\n";
	# parallel execution of BPP instances
	for (my $rep_i = 0; $rep_i < $repeats; $rep_i++) {
		my $wdir_rep;
		$wdir_rep = $wdir.'repeat'.$rep_i.'/';
		for (my $i = 0; $i < scalar @analyses; $i++) {
			sleep(1) while(scalar threads->list(threads::running) >= $maxThreads); # Limit threads running
			push(@threads, threads->create (\&runBPP, $wdir_rep, $analyses[$i]));  # Create and run next thread
			print scalar localtime() . ": Executed $analyses[$i] in repeat$rep_i\n";
			sleep(5); # wait some seconds - I think that strange things might happen otherwise...
		}
	}
	# Thread cleanup
	sleep(5);
	foreach my $thread (threads->list(threads::all)) {
		$thread->join();
	}
}



exit;
#=~ subs ~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
sub makeControlBPP {
	# get arguments
	my $wdir                = shift @_;
	my $analysis            = shift @_;
	my $SpeciesNumberString = shift @_;
	my $tree                = shift @_;
	my $thetaprior          = shift @_;
	my $tauprior            = shift @_;
	
	# make subdir
	print "Making subdirectory for analysis $analysis\n";
	mkdir ($wdir);
	chdir ($wdir);
	mkdir ($wdir.$analysis);
	chdir ($wdir.$analysis);

	# write file
	my $file = $analysis . '.ctl';
	print "Writing control file $file\n";
	open(my $FH, ">", $file) or die "Can't write to '$file': $!\n";
	
	print $FH "* BPP3 control file written by BPPmulti.pl\n";
	print $FH "
seed     = -1
	
seqfile  = ../../$seqfile
Imapfile = ../../$Imapfile
outfile  = $analysis.out.txt
mcmcfile = $analysis.mcmc.txt

speciesdelimitation = $speciesdelimitation
speciestree         = $speciestree
speciesmodelprior   = $speciesmodelprior

species&tree = $SpeciesNumberString
 $tree

usedata    = $useseqdata
nloci      = $nloci
cleandata  = $cleandata

thetaprior = $thetaprior
tauprior   = $tauprior

heredity   = $heredity ../../$heredityfile
locusrate  = $locusrate ../../$locusratefile

finetune   = $finetune

print      = $print
burnin     = $burnin
sampfreq   = $sampfreq
nsample    = $nsample
";
close $FH or die "Can't close '$file': $!\n";
}


sub makeControliBPP {
	# get arguments
	my $wdir                = shift @_;
	my $analysis            = shift @_;
	my $SpeciesNumberString = shift @_;
	my $tree                = shift @_;
	my $thetaprior          = shift @_;
	my $tauprior            = shift @_;
	
	# make subdir
	print "Making subdirectory for analysis $analysis\n";
	mkdir ($wdir);
	chdir ($wdir);
	mkdir ($wdir.$analysis);
	chdir ($wdir.$analysis);

	# write file
	my $file = $analysis . '.ctl';
	print "Writing control file $file\n";
	open(my $FH, ">", $file) or die "Can't write to '$file': $!\n";
	
	print $FH "* iBPP control file written by BPPmulti.pl\n";
	print $FH "
seed     = -1
	
seqfile  = ../../$seqfile
Imapfile = ../../$Imapfile
traitfile = ../../$traitfile
outfile  = $analysis.out.txt
mcmcfile = $analysis.mcmc.txt

speciesdelimitation = $speciesdelimitation
uniformrootedtrees  = $speciesmodelprior

species&tree = $SpeciesNumberString
 $tree

useseqdata   = $useseqdata
usetraitdata = $usetraitdata
nloci        = $nloci
ntraits      = $ntraits
nindT        = $nindT
cleandata    = $cleandata

thetaprior   = $thetaprior
tauprior     = $tauprior
nu0          = $nu0
kappa0       = $kappa0

heredity   = $heredity ../../$heredityfile
locusrate  = $locusrate ../../$locusratefile

finetune     = $finetune

print        = $print
burnin       = $burnin
sampfreq     = $sampfreq
nsample      = $nsample
";
close $FH or die "Can't close '$file': $!\n";
}

sub parseImap {
	# get arguments
	my $ImapFile = shift @_;
	
	my %number_of;
	
	open (my $FH, "<", $ImapFile) or die "Can't open '$ImapFile': $!\n";
	while (my $line = <$FH>) {
		$line =~ s/(^\s+)|(\s+$)//g;
		my @split = split(/\s+/, $line);
		$number_of{$split[1]}++;
	}
	
	my $SpeciesNumberString = keys %number_of;
	foreach my $key (sort (keys %number_of)) {
		$SpeciesNumberString = $SpeciesNumberString . ' ' . $key;
	}
	$SpeciesNumberString = $SpeciesNumberString . "\n";
	foreach my $key (sort (keys %number_of)) {
		$SpeciesNumberString = $SpeciesNumberString . ' ' . $number_of{$key};
	}
	
	return $SpeciesNumberString;
	close $FH or die "Can't close '$ImapFile': $!\n";
}

sub runBPP {
# the routine to start in the above multithread section:
	my $wdir     = shift @_;
	my $analysis = shift @_;
	chdir ($wdir.$analysis);
	
	# need to check wether this is an iBPP or an BPP analysis:
	my $cmd;
	if ($analysis =~ m/iBPP/) {
		$cmd = "$iBPP $wdir$analysis/$analysis.ctl > screen.$analysis.txt";
	} else {
		$cmd = "$BPP $wdir$analysis/$analysis.ctl > screen.$analysis.txt";
	}
	
	async{`$cmd`};
}
