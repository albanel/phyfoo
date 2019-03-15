#!/usr/bin/perl
use strict;
use warnings;
my $v = "salut toi";
substr($v,5,1) = "ation à ";
print("$v\n");
