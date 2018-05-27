/**
    PLCP
    Copyright (C) 2018 Lorraine A.K. Ayad and Carl Barton and Panagiotis  
    Charalampopoulos and Costas S. Iliopoulos and Solon P. Pissis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <limits.h>
#include <sys/time.h>
#include "plcp.h"


static struct option long_options[] =
 {
   { "input-file",              required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "hamming-dist",            required_argument, NULL, 'k' },
   { "threads",                 required_argument, NULL, 't' },
   { "rmqs",                	required_argument, NULL, 'r' },
   { "help",                    no_argument,       NULL, 'h' },
   {  NULL,                     0,                 NULL,  0  }
 };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> input_filename                 = NULL;
   sw -> output_filename                = NULL;
   sw -> k                              = 1;
   sw -> t                              = 1;
   sw -> r                              = 0;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "i:o:k:t:r:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
         case 'i':
           sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> input_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
           break;

         case 'k':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> k = val;
           args ++;
           break;

         case 't':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> t = val;
           break;

	 case 'r':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> r = val;
           break;
       }
    }

   if ( args < 3 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }


/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, " plcp <options>\n" );
   fprintf ( stdout, " Mandatory arguments:\n" );
   fprintf ( stdout, "  -i, --input-file	<str>	Sequence input filename (FASTA format).\n" );
   fprintf ( stdout, "  -o, --output-file	<str>	Output filename.\n" );
   fprintf ( stdout, "  -k, --hamming-dist	<int>	Hamming distance between matches.\n");
   fprintf ( stdout, " Optional arguments:\n" );
   fprintf ( stdout, "  -r, --rmqs		<int>	0 for computing rmqs using O(nlogn) space or 1 for O(n) space. (default: 0).\n");
   fprintf ( stdout, "  -t, --threads		<int>	The number of threads to be used (default: 1).\n\n");
 }

double gettime( void )
 {
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
 }

