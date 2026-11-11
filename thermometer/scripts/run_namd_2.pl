#!/usr/bin/perl
use warnings;
use strict;

# Run Minimization

my $folder = $ARGV[0];


opendir(my $dir, $folder) || die ("unbale to read $folder\n");
# deterministic order, and skip . / .. and dotfiles
my @files = sort grep { !/^\./ } readdir($dir);
close($dir);

foreach my $pdb (@files) {
    
    if ($pdb=~/.pdb$/) {
        
        my $psf= $pdb;
        $psf =~s/.pdb$/.psf/;
        print $psf . "\n";
        print $pdb . "\n";
        
        my $check = readConf($pdb, $psf);
        
        #my $command = "/home/mattia/Dropbox/Projects/THERMOSTABILITY/Mesofilization/ServerPipeline/NAMD_2.7b1_Linux-x86_64/charmrun /home/mattia/Dropbox/Projects/THERMOSTABILITY/Mesofilization/ServerPipeline/NAMD_2.7b1_Linux-x86_64/namd2 ++local +p2 /home/mattia/Dropbox/Projects/THERMOSTABILITY/Mesofilization/ServerPipeline/PDB-PSF/MIN_new.conf 2>> err.log";
	
	# Tool root comes from THERMO_HOME so NAMD-MIN can live in a per-job work
	# dir instead of two levels under the tool dir.
	my $home  = $ENV{THERMO_HOME} || "../..";
	my $procs = $ENV{NAMD_PROCS}  || 2;
	my $base  = $pdb; $base =~ s/\.pdb$//;

	my $command = "$home/bin/charmrun $home/bin/namd2 ++local +p$procs "
	            . "MIN_new.conf > '$base.namd.log' 2>&1";
	        
	print $command;
        # Do not die on failure: one bad structure must not abort the batch.
        # The caller detects it from the missing DCD / ENERGY csv.
        my $rc = system($command);
        warn "NAMD failed for $pdb (rc=$rc), see $base.namd.log\n" if $rc != 0;
        
    }
    
}


sub readConf {
    
    my $pdb = shift;
    my $psf = shift;
    
    my $pdbfolder= "pdb_dataset/";
    my $psffolder= "psf_dataset/";
    
    
    open(my $in, "<", "MIN_general.conf") || die ("dunable to read file\n");
    open(my $out, ">", "MIN_new.conf")|| die ("unable to open file\n");
    
    
    while (my $line=<$in>) {
        
        if ($line=~/^(structure\s+)/) {
            
            print "yooo!!!";
            print $out $1 . $psffolder . $psf . "\n";
            
            print $1 . $psf . "\n";
        }
        
        elsif ($line=~/^(coordinates\s+)/) {
            
            print $out $1 . $pdbfolder . $pdb . "\n";
            
            print $1 . $pdb . "\n";
            
        } #
        
        elsif ($line=~/^(outputName\s+)/) {
            
            # Was a 4-character truncation of the pdb name. That is THE root
            # cause of the broken chain: NAMD wrote e.g. "AF-O.dcd", catDCD.sh
            # then derived base="AF-O" and looked for "AF-O_autopsf.pdb", which
            # autopsf never wrote (it uses the full name). Use the full basename.
            # Capture $1 first: the substitutions below are themselves
            # successful matches and would reset it.
            my $prefix = $1;
            my $out_name = $pdb;
            $out_name =~ s/\.pdb$//;
            # psfgen feeds us <name>_autopsf.pdb; catDCD.sh re-derives the
            # structure file as <base>_autopsf.pdb, so the outputName base
            # must be <name>, not <name>_autopsf.
            $out_name =~ s/_autopsf$//;
            
            print $out $prefix . $out_name . "\n";
            
            print $prefix . $out_name . "\n";
            
        }
        
        else {
            
            print $out $line;
            
            
        }
        
    }
    
    close($in);
    close($out);
    
    return(1);
    
}
