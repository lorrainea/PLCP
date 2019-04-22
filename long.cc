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
#include <omp.h>
#include "plcp.h"

using namespace std;

INT long_plcp( unsigned char * x, struct TSwitch  sw, INT * PLCP, INT * P, INT * SA, INT * LCP, unordered_set< pair<INT,INT> , pair_hash > * h_map )
{	
	INT l = strlen( (char*) x);

	#pragma omp parallel
	for( unordered_set<pair<INT, INT>, pair_hash>::iterator it=h_map->begin(); it!=h_map->end(); ++it )
	{
		INT first = it->first;
		INT second = it->second;	
	
		compute_long( PLCP, P, sw, x, first, second );
	}
	
		
return 0;
}

INT compute_long( INT * PLCP, INT * P, TSwitch sw, unsigned char * x, INT first, INT second ) 
{

	INT l = strlen( (char*) x);
		    	
	INT match = 0;
	INT error = 0;

	while( first + match < l && second + match < l && error <= sw . k )
	{
		if( x[first + match] != x[second + match] )
		{
			error++	;
			
		}
		if( error <= sw . k )	
			match++;
	}

	if( match > PLCP[first] )	
	{
		PLCP[first] = match;
		P[first] = second;
	}
	if( match > PLCP[second] )
	{
		PLCP[second] =  match;
		P[second] = first;
	}

return 0;
}
