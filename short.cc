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
#include <algorithm>  
#include <vector>
#include <fstream>
#include <sys/time.h>

#include "plcp.h"
#include "short.h"

using namespace std;

unsigned int short_plcp( int i, char * alphabet, unsigned char * x, struct TSwitch  sw, unsigned int * PLCP, unsigned int * P, INT * SA, INT * LCP, INT * invSA )
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

	for(int c = 0; c<sw . k; c++)
		character[c] = 0;
	character[min( (int)sw . k , l - i )] = '\0';
	error_pos[min( (int) sw . k , l - i )] = '\0';

	int initial_error = sw .k ;
		
	if( initial_error > 0 )
		positions( alphabet, character, error_pos, i, initial_error, SA, PLCP, P, sw, x, start, end, len_ext);

	free( error_pos );
	free( character );
		
		
return 1;
}


unsigned int positions( char * alphabet, unsigned int * character, unsigned int * error_pos, int i, int error, INT *  SA, unsigned int * PLCP, unsigned int * P, struct TSwitch sw, unsigned char * x, int start, int end, int len_ext)
{ 
	for(int j=start; j<= end; j++) 
	{	
 		error_pos[sw . k - error ] = j;
		int length_k = 0;	

		for(int a = 0; a<strlen( (char*) alphabet ); a++)
		{
			if( alphabet[a] != x[j] )
			{
				character[sw . k - error] = a;

				int current_error =   sw . k - error;

				len_ext = extension( alphabet, character, error_pos, i, SA, PLCP, P, sw, x, current_error );	

				if( error  > 1 ) 
				{
					positions( alphabet, character, error_pos, i, error - 1, SA, PLCP, P, sw, x,  min( j+1, (int ) strlen( ( char * ) x )-1 ), min( j + len_ext+1, (int ) strlen( ( char * ) x )-1 ), len_ext);

				}
			}
		}

	}

return 1;
}

unsigned int extension( char * alphabet, unsigned int * character, unsigned int * error_pos, int i, INT * SA, unsigned int * PLCP, unsigned int * P, struct TSwitch  sw, unsigned char * x, int  error)
{

	int r  = strlen( (char*) x ) - 1;
	int l =  0;
	int len_ext = 0;
	int search = floor( ( r + l ) / 2.0 );

	bool end_reached = false;

	
	int R = sw . m + i + 5 ;

	while( r >= l  )
	{			
		
		int toAdd = 0;	
		int err = 0;
		int length = 0;

		for(int k=i; k< min (  R   , (int) strlen( (char*) x ) ); k++)
		{	
			char a = x[k];
			
			if( err <= error && error_pos[err] == k  )
			{
				a = alphabet[character[err]];
				err++;
			}

			if( a > x[SA[search]+toAdd] )
			{
				l = search + 1;
				search = floor( ( r + l ) / 2.0 );
				break;	
			}
			else if ( a < x[SA[search]+toAdd] )
			{
				r =  search - 1;
				search = floor( ( r + l ) / 2.0 );
				break;
			}
			else 
			{
				length++;
				if( length > len_ext )
				{
					len_ext = length;

					if( len_ext > PLCP[i] && SA[search] != i )
					{
						PLCP[i] =  len_ext;
						P[i] = SA[search];
					}
				}
			}

			if( k ==  min (  R - 1  , (int) strlen( (char*) x ) - 1 )  )
			{
				end_reached = true;
				break;
			}

			toAdd++;
		}

		if( end_reached == true )
			break;
	}
	
return len_ext;
}

