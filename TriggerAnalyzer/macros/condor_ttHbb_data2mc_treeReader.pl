#!/usr/local/bin/perl

####make DIRs
@dir = ("Condor", "Script", "Log", "Error", "Output");   

for $dir (@dir){
    if(! -d $dir ){
	mkdir $dir
    }
}


$script = "csv_normFix_treeReader";       # Name of job

$workingDir = "/uscms_data/d2/lwming/Trigger2016/CMSSW_8_0_8/src/TriggerRun2/TriggerAnalyzer";

$sample = 2500;
$intLumi = -1;
$tthf   = -1;
$Njobs  = 20;
$useHTbins = 0;

$num = @ARGV;

if( $num >= 1 ){
  $sample = $ARGV[0];

  if( $num>=2 ){
    $Njobs = $ARGV[1];

    if( $num>=3 ){
      $intLumi = $ARGV[2];
      
      if( $num>=4 ){
	$tthf = $ARGV[3];

	if( $num>=5 ){
	  $useHTbins = $ARGV[4];
	}
      }
    }
  }
}



open SHFILE, "> Script/condor\_$script.sh";
print SHFILE "#!/bin/sh\n";
print SHFILE "\n";
print SHFILE "echo $1\n";
print SHFILE "\n";
print SHFILE "echo \"\"\n";
print SHFILE "echo \"Using ROOT on Condor\"\n";
print SHFILE "echo \"\"\n";
print SHFILE "cd \${_CONDOR_SCRATCH_DIR}\n";
print SHFILE "\n";
print SHFILE "sample=\$1\n";
print SHFILE "NumEvents=\$2\n";
print SHFILE "NumJobs=\$3\n";
print SHFILE "jobN=\$4+1\n";
print SHFILE "intLumi=\$5\n";
print SHFILE "tthf=\$6\n";
print SHFILE "useHTbins=\$7\n";
print SHFILE "\n";
#print SHFILE "root -b -q $workingDir/macros/head.C '$workingDir/macros/$script.C+('\$sample','\$NumEvents','\$NumJobs','\$jobN','\$intLumi','\$tthf','\$useHTbins',1)'\n";
print SHFILE "root -b -q head.C $script.C'('\$sample','\$NumEvents','\$NumJobs','\$jobN','\$intLumi','\$tthf','\$useHTbins',1)'\n";
print SHFILE "\n";
close SHFILE;


open CONDORFILE, "> Condor/condor\_$script.jdl";
print CONDORFILE "# A Condor submission file\n";
print CONDORFILE "Executable              = Script/condor\_$script.sh\n";
print CONDORFILE "Universe                = vanilla\n";
print CONDORFILE "Getenv                  = true\n";
print CONDORFILE "\n";
print CONDORFILE "Arguments               = $sample -1 $Njobs \$(Process) $intLumi $tthf $useHTbins\n";
print CONDORFILE "Output                  = Output/condor\_$sample\_$script\_\$(Process)\_$tthf.out\n";
print CONDORFILE "Error                   = Error/condor\_$sample\_$script\_\$(Process)\_$tthf.err\n";
print CONDORFILE "Log                     = Log/condor\_$sample\_$script\_\$(Process)\_$tthf.log\n";
print CONDORFILE "\n";
print CONDORFILE "use_x509userproxy = true\n";
print CONDORFILE "Should_Transfer_Files   = YES\n";
print CONDORFILE "When_To_Transfer_Output = ON_EXIT\n";
print CONDORFILE "Transfer_Input_Files = $workingDir/data/csv_rwt_fit_hf_v2_final_2017_1_10test.root, $workingDir/data/csv_rwt_fit_lf_v2_final_2017_1_10test.root, $workingDir/data/cmva_rwt_fit_hf_v0_final_2017_1_10.root, $workingDir/data/cmva_rwt_fit_lf_v0_final_2017_1_10.root, $workingDir/macros/head.C, $workingDir/macros/$script.C \n";
print CONDORFILE "\n";
print CONDORFILE "#+IsLocalJob             = true\n";
print CONDORFILE "#Rank                    = TARGET.IsLocalSlot\n";
print CONDORFILE "\n";
print CONDORFILE "Queue $Njobs\n";
print CONDORFILE "\n";
close CONDORFILE;

system("chmod a+x Script/condor\_$script.sh");
print "submitting: condor_submit Condor/condor\_$script.jdl\n";
system("condor_submit Condor/condor\_$script.jdl");

