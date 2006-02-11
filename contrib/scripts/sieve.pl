#!/usr/bin/perl
#
##############################################################################
# sieve.pl                                                                   #
# Distributed sieving client script.                                         #
# Copyright 2006, Sten.                                                      #
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
#   1) Change $nick script parameter.  
#   2) Place these files into the same directory: 
#      - gnfs-lasieve4I14e.exe 
#      - sieve.pl 
#      - factorbase.dat (optional, just rename  block-XXX.job.afb.0 into factorbase.dat) 
#   3) Run the script. 
#   4) Ensure it works and forget about it for several weeks. :)
# 

use Socket;

sub ExecJob($);
sub GetAndParseJob($);

# --------------------------
# Configuration
# --------------------------
$nick = 'ChangeMe';             # !!!!!!!!!!!! change it !!!!!!!!!!!!!
$host = 'factoring.org';
$blocksize = 'small';

$content = "name=$nick&blocksize=$blocksize&blockcount=1";

$http = "POST /distributed/getblock.php HTTP/1.0\n".
        "Host: $host\n".
        "User-Agent: Internet Exploider 6.0\n".
        "Referer: $host\n".
        "Connection: close\n".
        "Content-Type: application/x-www-form-urlencoded\n".
        "Content-Length: ". length($content). "\n\n".
        $content;

$tcp = getprotobyname('tcp') or die "Couldn't getprotobyname!\n";
$hosti = inet_aton($host) or die "Couldn't look up host!\n";
$hosts = sockaddr_in(80, $hosti);

#
# Open .lastjob file and see if we have terminated unexpectedly last time
#
if (open(LASTJOBFILE, ".lastjob") && ($JobName = <LASTJOBFILE>))
{
    close (LASTJOBFILE);
    print "Restoring last job: $JobName\n";
    ExecJob($JobName);
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

    #
    # Call remote getblock.php script and allocate new task
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
    if (!($response =~ m|"http://$host(.*\.job)"|g))
    {
        print "Error: bad server response.\n";
        goto sleep_continue;
    }

    #
    # Download .job file and parse it, save to disk. 
    # N.B. Function sets $block, $nblocksize global vars as a side effect.
    #
    if (!($JobName = GetAndParseJob($1)))
    { 
        print "Error: GetAndParseJob() failed.\n";
        goto sleep_continue;
    }

    #
    # Execute .job file
    #
    ExecJob($JobName);

    #
    # Count number of relations we have found and add record to the log file
    #
    if (open(FILE, "< spairs-$block.out"))
    {
        $count = 0;
        $count++ while <FILE>;
        close (FILE);

        # save log line
        open (LOGFILE, '>>sieve.log');
        print LOGFILE "range:\t$block\tsize:\t$nblocksize\tnumrels:\t$count\n";
        close (LOGFILE);
    }
    else
    {
        print "Log Error!\n";
    }

    next;
 
sleep_continue:
    print "Error: Some error occured. Sleeping 10 seconds.\n";
    sleep(10); #wait 10 seconds and try again;
}

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
              "User-Agent: Internet Exploider 6.0\n".
              "Referer: $host\n".
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
    (($block) = ($response =~ /q0:\s+(\d+)/g)) || return 0; 
    (($nblocksize) = ($response =~ /qintsize:\s+(\d+)/g)) || return 0;

    #
    # Save job to file
    #
    open(JOBFILE, ">block-$block.job") || return 0; 
    print JOBFILE $response;
    close (JOBFILE);

    #
    # Save .lastjob file
    #
    open(LASTJOBFILE, ">.lastjob") || return 0; 
    print LASTJOBFILE "block-$block.job";
    close (LASTJOBFILE);
   
    return "block-$block.job";
}

sub ExecJob($)
{
    my $JobName = shift(@_);

    open(JOBFILE, $JobName) || die "Cann't open job file ($JobName).\n";

    @JobData = <JOBFILE>;
    close (JOBFILE);

    $JobData = join('', @JobData);

    #
    # Parse .job file and get q0 and qintsize parameters.
    #
    (($block) = ($JobData =~ /q0:\s+(\d+)/g)) || return 0; 
    (($nblocksize) = ($JobData =~ /qintsize:\s+(\d+)/g)) || return 0;
    print "block: $block, size: $nblocksize\n";

    #
    # Delete spairs-$block.out for current block (in case we were terminated
    # abnormally last time, now we restart last block from the beginning).
    #
    unlink "spairs-$block.out", ".last_spq1234";

    #
    # 1) Copy factorbase.dat to .afb.0 file. Dunno yet how to tell gnfs-lasieve4I14e
    #    to use given factorbase file.
    # 2) Execute lasieve
    # 3) Cleanup everything
    #
    if ($^O eq "MSWin32")
    {
        # Win32
        system ("copy factorbase.dat block-$block.job.afb.0");
        system ("gnfs-lasieve4I14e -k -o spairs-$block.out -v -a block-$block.job -n1234") == 0 
           || die "Error: Unable to execute gnfs-lasieve4I14e.\n";
    }
    else
    {
        # Unix
        system ("cp factorbase.dat block-$block.job.afb.0");
        system ("gnfs-lasieve4I14e -k -o spairs-$block.out -v -a block-$block.job -n1234") == 0 
           || die "Error: Unable to execute gnfs-lasieve4I14e.\n";
    }

    unlink ".lastjob", "block-$block.job", "block-$block.job.afb.0", ".last_spq1234", "ggnfs.log";
};
