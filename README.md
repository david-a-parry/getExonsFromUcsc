# getExonsFromUcsc

## Installation

To use this script you will need to use git to clone this repo and submodules. 

    git clone --recursive https://github.com/gantzgraf/getExonsFromUcsc.git

You will also need to install the following perl modules from CPAN if not 
already on your system:

    DBI
    DBD::mysql
    LWP::Simple

## Running

Run this script by either invoking it through 'perl getExonsFromUcsc.pl' or 
directly (e.g. ./getExonsFromUcsc.pl from the script directory). In its simplest form:

    ./getExonsFromUcsc.pl -g ABCD1

...will retrieve all exons from the gene ABCD1. 

To get details on usage and options, invoke with the --help or --manual options:
    
    ./getExonsFromUcsc.pl -h 

or 

    ./getExonsFromUcsc.pl --manual

The latter gives slightly more information and allows you to browse the 
documentation via a pager (press 'q' to exit).

### AUTHOR

David A. Parry

University of Edinburgh

### COPYRIGHT AND LICENSE

Copyright 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.


 
