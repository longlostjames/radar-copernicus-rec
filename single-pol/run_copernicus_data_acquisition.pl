#!/usr/bin/perl -w 
# Owain Davies 08062004

MAIN:
{

while(1)
{
	print "check tx power\n";
	$program = "/usr/local/bin/check_copernicus_power.pl";
	system($program);
	print "invoking recording program\n";
	$program = "/usr/local/bin/radar-copernicus-rec > /dev/null";
	system($program);
	print "sleeping for three seconds before restart\n";
	sleep 3;
}

}

