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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sys/time.h>

#include "plcp.h"

using namespace std;

unsigned int short_plcp( int i, char * alphabet, unsigned char * x, struct TSwitch  sw, unsigned int * PLCP, unsigned int * P, INT * SA, INT * LCP, INT * invSA, INT * A )
{	
	int l = strlen( (char*) x );

	int len = 0;

	if( invSA[i]+1 < l )
		len = max( LCP[invSA[i]], LCP[invSA[i]+1] ) + 1;
	else len  = LCP[invSA[i]] + 1;

	int len_ext = len;

	int start = i;
	int end = min(i + len, l-1 );
	unsigned int * error_pos = ( unsigned int * ) calloc( min( (int) sw . k + 1, l - i + 1 )  , sizeof( unsigned int ) );
	unsigned int * character = ( unsigned int * ) calloc( min( (int) sw . k + 1, l - i + 1 )  , sizeof( unsigned int ) );

	unsigned int * prev_errors = ( unsigned int * ) calloc( min( (int) sw . k + 1, l - i + 1 )  , sizeof( unsigned int ) );

	for(int c = 0; c<sw . k; c++)
		character[c] = 0;
	character[min( (int)sw . k , l - i )] = '\0';
	error_pos[min( (int) sw . k , l - i )] = '\0';

	int initial_error = sw .k;
		
	if( initial_error > 0 )
		positions( alphabet, character, error_pos, i, initial_error, SA, invSA, LCP, PLCP, P, sw, x, start, end, len_ext, prev_errors, A);

	free( error_pos );
	free( character );
	free( prev_errors );
		
	return 1;
}


unsigned int positions( char * alphabet, unsigned int * character, unsigned int * error_pos, int i, int error, INT *  SA, INT * invSA, INT * LCP, unsigned int * PLCP, unsigned int * P, struct TSwitch sw, unsigned char * x, int start, int end, int len_ext,  unsigned int * prev_errors, INT * A)
{ 
	for(int j = start; j <= end; j++) 
	{	
 		error_pos[sw . k - error ] = j;
		
		int length_k = 0;	

		for(int a = 0; a<strlen( (char*) alphabet ); a++)
		{
			if( alphabet[a] != x[j] )
			{
				character[sw . k - error] = a;

				int current_error =   sw . k - error;

				len_ext = extension( alphabet, character, error_pos, i, SA, invSA, LCP, PLCP, P, sw, x, current_error, prev_errors, A  );	

				if( error  > 1 ) 
				{
					positions( alphabet, character, error_pos, i, error - 1, SA, invSA,  LCP, PLCP, P, sw, x,  min( j+1, (int ) strlen( ( char * ) x )-1 ), min( j + len_ext+1, (int ) strlen( ( char * ) x )-1 ), len_ext, prev_errors, A );

				}
			}
		}

	}

	return 1;
}

void binary_search ( INT * left, INT * right, INT * key, INT value )
{
	if ( value )
	{
		( * left ) = ( * key ) + 1;
		( * key ) = floor( ( (* right ) + ( * left ) ) / 2.0 );
	}
	else
	{
		( * right )  =  ( * key ) - 1;
		( * key ) = floor( ( ( * right ) + ( * left ) ) / 2.0 );
	}
}

unsigned int extension( char * alphabet, unsigned int * character, unsigned int * error_pos, int i, INT * SA, INT * invSA, INT * LCP, unsigned int * PLCP, unsigned int * P, struct TSwitch  sw, unsigned char * x, int  error, unsigned int * prev_errors,  INT * A)
{
	INT error_p = error_pos[ error ];
	INT prev_error = 0;

	if( error == 0 )		prev_error = i;
	else				prev_error = prev_errors[ error - 1];

	INT cur_err = prev_error + error_p - i;

	INT n  = strlen( (char*) x );
	INT r  = n - 1;

	INT l =  0;
        INT len_ext = 0;
	INT search = floor( ( r + l ) / 2.0 );
	while( r >= l )
	{
	        INT rq = 0;
		INT rq2 = 0;
		int pos_value = SA[search];
		 
		if ( invSA[prev_error] == search ) rq = n - prev_error;
		else rq = LCP[rmq(A, LCP, n , invSA[prev_error], search )];
//cout<<SA[search]<<"  ";
		if( rq >= error_p - i ) // r - j 
		{	
			if( SA[search] + error_p - i == n ) binary_search ( &l, &r, &search, 1 );	//q reached end of string
			else if( alphabet[character[error]] == x[SA[search] + error_p - i] )	//letter comparison successful
			{
				if( error_p<n-1 && SA[search] + error_p - i <n-1)	//RMQ allowed, else we are good, rq2=0
				{
					if (invSA[error_p +1] == invSA[SA[search] + error_p - i +1] ) rq2 = n - ( error_p +1 );
					else rq2 = LCP[rmq(A, LCP, n , invSA[error_p +1], invSA[SA[search] + error_p - i +1] )];
				}

				INT value = 0; if( SA[search] +  error_p - i + rq2 == n-1 || ( error_p + rq2 < n-1 && x[ error_p + rq2 + 1 ] > x[ SA[search] +  error_p - i + rq2 + 1 ]) ) value = 1;
				binary_search ( &l, &r, &search, value );
	
				if(len_ext<error_p - i + rq2 + 1) 
				{
					len_ext = error_p - i+ rq2 + 1;
					prev_errors[error]=pos_value;
					if( len_ext > PLCP[i] && pos_value != i )
					{
						PLCP[i] = len_ext;
						P[i] = pos_value;
					}		
				}
			}
			else //letter comparison unsuccessful
			{
			        INT value = 0; if (  alphabet[character[error]] > x[SA[search]+error_p-i]  ) value = 1;
				binary_search ( &l, &r, &search, value );
			}
		}
		else //( rq < error_p - i )
		{	
			INT value = 0; if ( SA[search] + rq == n || x[ prev_error + rq ] > x[SA[search] + rq]  ) value = 1;
		        binary_search ( &l, &r, &search, value  );
		}
//cout<<rq<<"  "<<alphabet[character[error]]<<"  "<<error_p<<"  "<<rq2<<endl;
	}//end of while

	return ( len_ext - error_p + i - 1 );
}

