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
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <string>
#include <sys/time.h>
#include <math.h>
#include <omp.h>
#include "plcp.h"

using namespace std;

int main( int argc, char **argv )
{

	struct TSwitch  sw;

	FILE *           in_fd;                  // the input file descriptor
	FILE *           out_fd;                 // the input file descriptor
        char *           input_filename;         // the input file name
        char *           output_filename;        // the output file name 
        unsigned char * seq       = NULL;          // the sequence in memory
        unsigned char * seq_id    = NULL;          // the sequence in memory

	/* Decodes the arguments */
        INT i = decode_switches ( argc, argv, &sw );

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
	}


	omp_set_num_threads( sw.t );

	double start = gettime();

        /* Read the file in memory */
        if ( ! ( in_fd = fopen ( input_filename, "r") ) )
        {
                fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
                return ( 1 );
        }

        INT max_alloc_seq_id = 0;
        INT max_alloc_seq = 0;
	INT seq_len = 0;

        char c;
	c = fgetc( in_fd );
        do
        {
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
			return ( 1 );
		}
		else
		{
			INT max_alloc_seq_id = 0;
			INT seq_id_len = 0;
			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id )
				{
					seq_id = ( unsigned char * ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id += ALLOC_SIZE;
				}
				seq_id[ seq_id_len++ ] = c;
			}
			seq_id[ seq_id_len ] = '\0';
			
		}
                
                INT max_alloc_seq_len = 0;

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

			seq[ seq_len++ ] = c;
                }

		seq[ seq_len ] = '\0';

        } 
	while( c != EOF );

	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}
	fprintf( stderr, "Processing sequence %s\n", seq_id);

	INT l = seq_len;

	INT * SA = ( INT * ) malloc( ( l ) * sizeof( INT ) );
	INT * LCP = ( INT * ) calloc  ( l, sizeof( INT ) );
	INT * invSA = ( INT * ) calloc( l , sizeof( INT ) );

	compute_SA( seq, l, SA );
	compute_invSA( seq, l, SA, invSA );
	compute_LCP( seq, l, SA, invSA, LCP );
	INT * A;

	if( sw . r == 0 )
	{
		INT lgl = flog2( l );
		A = ( INT * ) calloc( ( INT ) l * lgl, sizeof(INT) );
		rmq_preprocess(A, LCP, l);
	}

	INT * PLCP = ( INT * ) calloc( ( seq_len + 1 ) , sizeof( INT ) );
	PLCP[ seq_len ] = '\0';

	INT * P = ( INT * ) calloc( ( seq_len + 1 ) , sizeof( INT ) );
	P[ seq_len ] = '\0';

	populate_PLCP( seq, l, SA, invSA, LCP, PLCP, P );

	sw . m =  ( 2.0 ) *  log( l );

	unordered_map< pair<INT, INT>, INT,pair_hash > h_map;

	if ( sw . k >= l )
	{
		fprintf( stderr, " Error: k is too large!\n");
		return ( 1 );
	}

	if( sw . k > 0 )
	{
		fprintf( stderr, " Checking for short LCPs with %ld errors.\n", sw . k);
		short_plcp( seq, sw, PLCP, P, SA, LCP, invSA, A, &h_map );

		fprintf( stderr, " Checking for long LCPs with %ld errors.\n", sw . k);
		long_plcp( seq, sw, PLCP, P, SA, LCP, &h_map );
	}

	double end = gettime();


	/*for(INT i =0; i<l; i++)
		cout<<PLCP[i]<<" ";

	cout<<endl;

	for(INT i =0; i<l; i++)
		cout<<P[i]<<" ";*/

	if ( ! ( out_fd = fopen ( output_filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
		return ( 1 );
	}


	for( INT i=0; i < l; i++ ) 
		fprintf ( out_fd, "%ld ", PLCP[i] );

	fprintf( out_fd, "\n" );

	for( INT i=0; i < l; i++ ) 
		fprintf ( out_fd, "%ld ", P[i] );

	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

        fprintf( stderr, "Elapsed time: %lf secs\n", ( end - start ) );

	/* De-allocate */
	free ( PLCP );
	free ( P );
        free ( seq );
        free ( seq_id );
        free ( sw . input_filename );
        free ( sw . output_filename );
	free ( invSA );
	free ( SA );
	free ( LCP );

	return ( 0 );
}
