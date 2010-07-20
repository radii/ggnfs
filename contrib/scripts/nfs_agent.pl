#!/usr/bin/perl
#
##############################################################################
# sieve.pl                                                                   #
# Distributed sieving client script.                                         #
# Copyright 2006-2010, Sten.                                                 #
##############################################################################
#
#  This file is part of GGNFS.
#
#   GGNFS is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   GGNFS is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with GGNFS; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# How to use:
#
#   1) Change $nick script parameter and optionally set $output_dir to some shared folder (*.spairs output files will be moved there).  
#   2) Place these files into the same directory: 
#      - gnfs-lasieve4I14e.exe 
#      - sieve.pl 
#      - factorbase.dat (optional, just rename  block-XXX.job.afb.0 into factorbase.dat) 
#        !!! commands to copy factobase.dat tto block-XXX.job.afb.0 are currently disabled. !!!
#   3) Run the script. 
#   4) Ensure it works and forget it for several weeks. :)
#
use Socket;

sub ExecJob($);
sub GetAndParseJob($);
sub MoveSpairsToOutputDir();

# --------------------------
# Main configuration
# --------------------------
$nick       = 'ChangeMe';              # ex. 'Dmit';
$host       = 'factoring.org';         # http header field
$script     = 'sieve/getblock.php';    # ex. '2p643m3/getblock.php';
$task       = 'TaskName.1';            # Task name
$blocksize  = 'small';
$output_dir = undef;                   # ex. for Windows: "\\\\myserver\\spairs"; Output files will be automatically moved to \\myserver\spairs.

# --------------------------
# Proxy configuration
# --------------------------
$server    = 'factoring.org';          # ex. 'my.proxy.com';
$port      = 80;                       # ex. 3128
$hdr_add   = '';                       # ex. 'Proxy-Authorization: Basic a384ZrRab5Z7O9luchJKZdZ5ejB2QQ==\n';

$content = "name=$nick&blocksize=$blocksize&blockcount=1&task=$task";

$http = "POST http://$host/$script HTTP/1.0\n".
        "Host: $host\n".
        "User-Agent: Factoring Machine 1.0\n".
        "Referer: $host\n".
        $hdr_add.
        "Connection: close\n".
        "Content-Type: application/x-www-form-urlencoded\n".
        "Content-Length: ". length($content). "\n\n".
        $content;

$tcp = getprotobyname('tcp') or die "Couldn't getprotobyname!\n";
$hosti = inet_aton($server) or die "Couldn't look up host!\n";
$hosts = sockaddr_in($port, $hosti);

MoveSpairsToOutputDir();

#
# Open .lastjob file and see if we have terminated unexpectedly last time
#
if (open(LASTJOBFILE, ".lastjob") && ($jobFileName = <LASTJOBFILE>))
{
    close (LASTJOBFILE);
    print "Restoring last job: $jobFileName\n";
    ExecJob($jobFileName);
}

while (1) 
{
    #
    # check if we were signalled to stop processing
    #  
    if (-e '.stop')
    {
        unlink '.stop';
        last;
    }

    MoveSpairsToOutputDir();

    # 
    # Check TODO list 
    # 
    if (opendir(DIR, "todo")) 
    {
        while ($entry = readdir( DIR ))
        {
           if (not($entry eq '.' || $entry eq '..')) {last;}
        }

        if ($entry) 
        {
            $oldname = "todo\\$entry"; 
            $newname = $entry; 

            print "Renaming $oldname to $newname\n"; 
            rename ($oldname, $newname); 
            ExecJob ($newname); 
            next;
        }
    }

    #
    # Call remote getblock.php script and allocate new job
    #
    socket(SOK, PF_INET, SOCK_STREAM, $tcp)
        || goto sleep_continue;
    connect(SOK, $hosts)
        || goto sleep_continue;

    select SOK; $| = 1; select STDOUT;
    print SOK $http;

    $response = '';
    while (<SOK>) { $response .= $_; }
    close SOK;

    #
    # Extract .job file URL
    #
    if (!($response =~ m|"(http://$host.*\.job)"|g))
    {
        print "Error: bad server response.\n";
        goto sleep_continue;
    }

    #
    # Download .job file and parse it, save to disk. 
    # N.B. Function sets $blkStart, $blkSize, $blkEnd global 
    #      vars as a side effect.
    #
    if (!($jobFileName = GetAndParseJob($1)))
    { 
        print "Error: GetAndParseJob() failed.\n";
        goto sleep_continue;
    }

    #
    # Execute .job file
    #
    ExecJob($jobFileName);

    next;
 
sleep_continue:
    print "Error: Some error occured. Sleeping 10 seconds.\n";
    sleep(10); #wait 10 seconds and try again;
}

############################################################################
############################################################################
############################################################################

sub GetAndParseJob($)
{
    my $Url = shift(@_);

    #
    # Retrive .job file from the server
    #
    socket(SOK, PF_INET, SOCK_STREAM, $tcp)
        || return 0;
    connect(SOK, $hosts)
        || return 0;

    select SOK; $| = 1; select STDOUT;
    print SOK "GET $Url HTTP/1.0\n".
              "Host: $host\n".
              "User-Agent: Factoring Machine 6.0\n".
              "Referer: $host\n".
              $hdr_add. 
              "Connection: close\n".
              "Content-Type: application/x-www-form-urlencoded\n".
              "Content-Length: 0\n\n";

    $response = '';
    while (<SOK>) { $response .= $_; }
    close SOK;

    #
    # Remove HTTP header
    # 
    (($response) = ($response =~ /HTTP.*200.*\r\n\r\n(.*)/s)) || return 0;

    #
    # Parse .job file and get q0 and qintsize parameters. Also serves as
    # rudimentary job file syntax check.
    #
    (($baseName) = ($response =~ /name:\s+(\S+)/g)) || return 0;
    (($blkStart) = ($response =~ /q0:\s+(\d+)/g)) || return 0; 
    (($blkSize) = ($response =~ /qintsize:\s+(\d+)/g)) || return 0;
#    $blkEnd = $blkStart + $blkSize;

    #
    # Save job to file
    #
    $jobFileName = sprintf ("%s.%09d.%d.sieve.job", $baseName, $blkStart, $blkSize);
    open(JOBFILE, ">$jobFileName") || return 0; 
    print JOBFILE $response;
    close (JOBFILE);

    #
    # Save .lastjob file
    #
    open(LASTJOBFILE, ">.lastjob") || return 0; 
    print LASTJOBFILE "$jobFileName";
    close (LASTJOBFILE);
   
    return "$jobFileName";
}

############################################################################
############################################################################
############################################################################

sub ExecJob($)
{
    my $jobFileName = shift(@_);

    open(JOBFILE, $jobFileName) || die "Can't open job file ($jobFileName).\n";

    @JobData = <JOBFILE>;
    close (JOBFILE);

    $JobData = join('', @JobData);

    #
    # Parse .job file and get q0 and qintsize parameters.
    #
    (($baseName) = ($JobData =~ /name:\s+(\S+)/g)) || return 0;
    (($blkStart) = ($JobData =~ /q0:\s+(\d+)/g)) || return 0;
    (($blkSize) = ($JobData =~ /qintsize:\s+(\d+)/g)) || return 0;
#    $blkEnd = $blkStart + $blkSize;
    print sprintf ("block: %09d, size: %d\n", $blkStart, $blkSize);

    #
    # Delete $baseName.$blkStart-$blkEnd.work for current block (in case we were terminated
    # abnormally last time, now we restart last block from the beginning).
    #
    $workFileName = sprintf ("%s.%09d.%d.work", $baseName, $blkStart, $blkSize);
    $spairsFileName = sprintf ("%s.%09d.%d.spairs", $baseName, $blkStart, $blkSize);
    $spqFileName = ".last_spq1234";
    unlink "$workFileName", "$spqFileName";

    #
    # 1) Copy factorbase.dat to .afb.0 file. Dunno yet how to tell gnfs-lasieve4I14e
    #    to use given factorbase file.
    # 2) Execute lasieve
    # 3) Cleanup everything
    #
    if ($^O eq "MSWin32") 
    {
        # Win32
        system ("copy factorbase.dat $jobFileName.afb.0");
    }
    else {
        # Unix
        #system ("cp factorbase.dat $jobFileName.afb.0");
    }

    system ("gnfs-lasieve4I14e -k -o $workFileName -v -a $jobFileName -n1234") == 0 
       || die "Error: Unable to execute gnfs-lasieve4I14e.\n";

    rename ($workFileName, $spairsFileName);
    
    unlink ".lastjob", "$jobFileName", "$jobFileName.afb.0", "$spqFileName", "ggnfs.log";

    #
    # Count number of relations we have found and add record to the log file
    #
    if (open(FILE, "< $spairsFileName")) 
    {
        $count = 0;
        $count++ while <FILE>;
        close (FILE);

        # save log line
        open (LOGFILE, ">> $baseName.sieve.log");
        print LOGFILE sprintf ("block: %09d\tsize: %d\tnumrels: %d\n", $blkStart, $blkSize, $count);
        close (LOGFILE);
    }
    else
    {
        print "Log Error!\n";
    }
};

sub MoveSpairsToOutputDir()
{
    # 
    # Move .spairs files to the output directory.
    # 
    if (defined($output_dir))
    {
        if (opendir(DIR, "."))
        {       
            while (my $fileName = readdir(DIR)) 
            { 
                if ($fileName eq '.' || $fileName eq '..')
                {
                    next;
                }

                if (!($fileName =~ m/.*\.spairs/))
                {
                    next;
                }

                if ($^O eq "MSWin32") 
                {
                    # Win32
                    system ("move $fileName \"$output_dir\"");
                }
                else 
                {
                    # Unix
                    system ("mv $fileName \"$output_dir\"");
                }
            }

            closedir(DIR);
        }                
    }
}
