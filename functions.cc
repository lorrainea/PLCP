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
#include <fstream>
#include <sys/time.h>
#include "plcp.h"

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                           // include header for suffix sort
#endif


#include <sdsl/bit_vectors.hpp>	

using namespace std;

unsigned int LCParray( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP )
{										
	INT i=0, j=0;

	LCP[0] = 0;
	for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
		if ( ISA[i] != 0 ) 
		{
			if ( i == 0) j = 0;
			else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
			while ( text[i+j] == text[SA[ISA[i]-1]+j] )
				j++;
			LCP[ISA[i]] = j;
		}

return 1;
}


unsigned int compute_SA( unsigned char *text, INT n, INT * SA )
{

        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }

	#ifdef _USE_64
        if( divsufsort64( text, SA, n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

	#ifdef _USE_32
        if( divsufsort( text, SA,  n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

return 1;
}


unsigned int compute_invSA( unsigned char *text, INT n, INT * SA, INT * invSA )
{
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < n; i ++ )
        {
                invSA [SA[i]] = i;
        }

return 1;
}


unsigned int compute_LCP( unsigned char *text, INT n, INT * SA, INT * invSA, INT * LCP )
{
        if( ( LCP == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        if( LCParray( text, n, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }

return 1;
}



unsigned int populate_PLCP( unsigned char *text, int l, INT * SA, INT * invSA, INT * LCP, unsigned int * PLCP, unsigned int * P )
{

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

return 1;
}

