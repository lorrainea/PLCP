/**
    PLCP
    Copyright (C) 2017 Lorraine A.K. Ayad and Panagiotis Charalampopoulos 
    and Costas S. Iliopoulos and Solon P. Pissis

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
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <string>
#include <sys/time.h>
#include <math.h>
#include "short.h"
#include "plcp.h"
#include <omp.h>

using namespace std;

int main( int argc, char **argv )
{

	struct TSwitch  sw;

	FILE *           in_fd;                  // the input file descriptor
	FILE *           out_fd;                 // the input file descriptor
        char *           input_filename;         // the input file name
        char *           output_filename;        // the output file name
	char *          alphabet;  
	int i;
        unsigned char * seq    = NULL;          // the sequence in memory

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	/* Check the arguments */
        if ( i < 1 )
        {
                usage ();
                return ( 1 );
        }
        else
        {
                input_filename          = sw . input_filename;
                output_filename         = sw . output_filename;

		alphabet = ( char * ) DNA;
               /* else if ( ! strcmp ( "PROT", sw . alphabet ) )  alphabet = ( char * ) PROT;
                else
                {
                        fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences or `USR' for sequences over a user-defined alphabet!\n" );
                        return ( 1 );
                }*/
        }

	omp_set_num_threads( sw.T );

	double start = gettime();

        /* Read the file in memory */
        if ( ! ( in_fd = fopen ( input_filename, "r") ) )
        {
                fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
                return ( 1 );
        }

        char c;
        unsigned int num_seqs = 0;           
        unsigned int total_length = 0;       
        unsigned int max_alloc_seq_id = 0;
        unsigned int max_alloc_seq = 0;
	unsigned int seq_len = 0;

        do
        {
                
		if ( num_seqs >= max_alloc_seq )
                {
                        seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char  ) );
                        max_alloc_seq += ALLOC_SIZE;
                }

   
                unsigned int max_alloc_seq_len = 0;


                while ( ( c = fgetc( in_fd ) ) != EOF )
                {
                        if( seq_len == 0 && c == '\n' )
                        {
                                fprintf ( stderr, " Omitting empty sequence in file %s!\n", input_filename );
                                c = fgetc( in_fd );
                                break;
                        }
                        if( c == '\n' || c == ' ' ) 
				continue;
			
                        c = toupper( c );

                        if ( seq_len >= max_alloc_seq_len )
                        {
                                seq= ( unsigned char * ) realloc ( seq,   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
                                max_alloc_seq_len += ALLOC_SIZE;
                        }

			if( strchr ( alphabet, c ) )
			{
				seq[ seq_len++ ] = c;
			}
			else
			{
				fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", input_filename, c );
				return ( 1 );
			
                        }

                }

		seq[ seq_len ] = '\0';
		total_length += seq_len;
		num_seqs++;

        } 
	while( c != EOF );

	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	INT l = seq_len;

	INT * SA = ( INT * ) malloc( ( l ) * sizeof( INT ) );
	INT * LCP = ( INT * ) calloc  ( l, sizeof( INT ) );
	INT * invSA = ( INT * ) calloc( l , sizeof( INT ) );
	

	compute_SA( seq, l, SA );
	compute_invSA( seq, l, SA, invSA );
	compute_LCP( seq, l, SA, invSA, LCP );


	unsigned int * PLCP = ( unsigned int * ) calloc( ( seq_len + 1 ) , sizeof( unsigned int ) );
	PLCP[ seq_len ] = '\0';

	unsigned int * P = ( unsigned int * ) calloc( ( seq_len + 1 ) , sizeof( unsigned int ) );
	P[ seq_len ] = '\0';

	populate_PLCP( seq, l, SA, invSA, LCP, PLCP, P );

	int alph_len = 4;
	
	sw . m =  ceil( (sw.k+2) * ( log(l)/ log(alph_len) ) );

	if ( sw . k >= l )
	{
		fprintf( stderr, " Error: k is too large!\n");
		return ( 1 );
	}

	k_mappability( seq, sw, PLCP, P, SA, LCP  );

	#pragma omp parallel for
	for(int i =0; i<l; i++)
	{
		if( PLCP[i] < sw . m )
			short_plcp( i, alphabet, seq, sw, PLCP, P, SA, LCP, invSA);		
	}

	double end = gettime();

	if ( ! ( out_fd = fopen ( output_filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
		return ( 1 );
	}


	for(int i=0; i<seq_len; i++ ) fprintf ( out_fd, "%d ", PLCP[i] );

	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

        fprintf( stderr, "Elapsed time for comparing sequences: %lf secs\n", ( end - start ) );

	/* De-allocate */
	free ( PLCP );
	free ( P );
        free ( seq );
        free ( sw . input_filename );
        free ( sw . output_filename );
	free ( invSA );
	free ( SA );
	free ( LCP );

	return ( 0 );
}
