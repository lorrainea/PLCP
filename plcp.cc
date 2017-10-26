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
#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                           // include header for suffix sort
#endif

#include <sdsl/bit_vectors.hpp>	
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
	
	

	//cout<<x<<endl;
	INT * SA;
	INT * LCP;
	INT * invSA;
	INT l = seq_len;

        /* Compute the suffix plcp */
        SA = ( INT * ) malloc( ( l ) * sizeof( INT ) );
        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }

	#ifdef _USE_64
        if( divsufsort64( seq, SA, l ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

	#ifdef _USE_32
        if( divsufsort( seq, SA,  l ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

        /* Compute the inverse SA plcp */
        invSA = ( INT * ) calloc( l , sizeof( INT ) );
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < l; i ++ )
        {
                invSA [SA[i]] = i;
        }

	LCP = ( INT * ) calloc  ( l, sizeof( INT ) );
        if( ( LCP == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        /* Compute the LCP plcp */
        if( LCParray( seq, l, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }

	//cout<<"SA"<<endl;
	unsigned int * PLCP = ( unsigned int * ) calloc( ( seq_len + 1 ) , sizeof( unsigned int ) );
	PLCP[ seq_len ] = '\0';

	unsigned int * P = ( unsigned int * ) calloc( ( seq_len + 1 ) , sizeof( unsigned int ) );
	P[ seq_len ] = '\0';

	for(int i =0; i<l; i++)
	{
		if( invSA[i]+1 < l )
		{
			PLCP[i] = max( LCP[invSA[i]], LCP[invSA[i]+1]);
	
			if( PLCP[i] == LCP[invSA[i]]  )
				P[i] = SA[invSA[i]-1];
			else P[i] = SA[invSA[i]+1];
		}
		else 
		{
			PLCP[i] =  LCP[invSA[i]];
			P[i] = SA[invSA[i]-1];
		}
	}

	int alph_len = 5;
	
	sw . m =  ceil( (sw.k+2) * ( log(l)/ log(alph_len) ) );

	if ( sw . k >= l )
	{
		fprintf( stderr, " Error: k is too large!\n");
		return ( 1 );
	}
	for(int i =0; i<l; i++)
	{
		if( PLCP[i] >= sw . m )
			k_mappability( i, seq, sw, PLCP, P, SA, LCP  );
		else 
		{
			short_plcp( i, alphabet, seq, sw, PLCP, P, SA, LCP, invSA);
		}		
	}

	double end = gettime();

	/*if ( ! ( out_fd = fopen ( output_filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
		return ( 1 );
	}*/

	cout<<" PLCP: ";
	for(int i=0; i<seq_len; i++ )
		cout<<PLCP[i]<<" ";

		cout<<endl;
	
	cout<<" P: ";
	for(int i=0; i<seq_len; i++ )
		cout<<P[i]<<" ";
	cout<<endl;

	/*if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}*/

        fprintf( stderr, "Elapsed time for comparing sequences: %lf secs\n", ( end - start ) );

	/* De-allocate */
	free ( PLCP );
	free ( P );
        free ( seq );
        free ( sw . input_filename );
        //free ( sw . output_filename );
	free ( invSA );
	free ( SA );
	free ( LCP );

	return ( 0 );
}
