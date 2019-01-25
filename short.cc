/**
    PLCP
    Copyright (C) 2018 Lorraine A.K. Ayad and Panagiotis Charalampopoulos 
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
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <omp.h>
#include <sys/time.h>

#include "plcp.h"

using namespace std;

INT nchoosek( INT n, INT k ) 
{
	if (n < k) 
		return 0;
        if (n == k) 
		return 1;
        INT s = min(k, (n - k));
        INT t = n;
        for (INT a = n - 1, b = 2; b <= s; a--, b++)
            t = (t * a) / b;

return t;
}

INT combinadic( INT n, INT k, INT i, INT ** error_pos ) 
{
        INT x = i;
        INT a = n;
        INT b = k;

	INT thread = omp_get_thread_num();

        for (INT c = 0; c < k; c++) 
	{
		INT v = a - 1;
        	while ( nchoosek( v, b) > x ) 
			v--;
		
		error_pos[thread][c] = v;
           	x = x - nchoosek(error_pos[thread][c], b);
           	a = error_pos[thread][c];
           	b--;
	}

	for (INT j = 0; j < k; j++)
            	error_pos[thread][j] = n - error_pos[thread][j] - 1;

return 1;
}


INT short_plcp( unsigned char * x, TSwitch sw, INT * PLCP, INT * P, INT * SA, INT * LCP, INT * invSA, INT * A, unordered_map<pair<INT, INT>, INT, pair_hash> * h_map )
{	
	INT l = strlen( (char*) x );
	INT nck = nchoosek( sw.m, sw.k );
	
	INT ** thread_plcp = ( INT ** ) calloc( ( sw . t ) , sizeof( INT * ) ); 
	INT ** thread_p = ( INT ** ) calloc( ( sw . t ) , sizeof( INT * ) ); 
	INT ** column_lengths = ( INT ** ) calloc( ( sw.t) , sizeof( INT *) );
	INT ** start_column = ( INT ** ) calloc( ( sw.t ) , sizeof( INT *) );
	INT ** order = ( INT ** ) calloc( ( sw.t ) , sizeof( INT *) );
	INT ** rank1 = ( INT ** ) calloc( ( sw . t ) , sizeof( INT *) ); 
	INT ** rank2 = ( INT ** ) calloc( ( sw . t ) , sizeof( INT *) ); 
	INT ** final_rank = ( INT ** ) calloc( ( sw . t ) , sizeof( INT *) ); 
	INT ** error_pos = ( INT ** ) calloc( ( sw . t ) , sizeof( INT * ) ); 

	for(INT i =0; i<sw . t; i++)
	{
		thread_plcp[i] = ( INT * ) calloc( ( l ) , sizeof( INT  ) ); 
		thread_p[i] = ( INT * ) calloc( ( l ) , sizeof( INT ) );
		order[i] = ( INT * ) calloc( ( l ) , sizeof( INT ) );
		rank1[i] = ( INT * ) calloc( ( l ) , sizeof( INT ) );
		rank2[i] = ( INT * ) calloc( ( l ) , sizeof( INT ) );
		final_rank[i] = ( INT * ) calloc( ( l ) , sizeof( INT ) );
		column_lengths[i] = ( INT * ) calloc( ( sw.k + 1 ) , sizeof( INT ) );
		start_column[i] = ( INT * ) calloc( ( sw.k + 1 ) , sizeof( INT ) );
		error_pos[i] = ( INT * ) calloc( ( sw.k ) , sizeof( INT ) );	
	}

	vector<vector<vector<INT>>> * bucket = new vector<vector<vector<INT>>>;
	vector<vector<INT>> inp;
	vector<INT> inp2;

	for(INT a = 0; a<=l; a++)
		inp.push_back( inp2 ); // <=l as rank starts from 1 and 0 means no substring at all in column
	
	for( INT b = 0; b<sw.t; b++)
		bucket->push_back( inp );

	INT thread = omp_get_thread_num();


	#pragma omp parallel for	
	for( INT i =0; i<nck; i++)
	{
		combinadic( sw.m, sw.k, i, error_pos );
		compute_plcp( l, error_pos, SA, invSA, LCP, PLCP, P, A, thread_plcp, thread_p, sw, column_lengths, start_column, order, rank1, rank2, final_rank, bucket, h_map );
	}


	for( INT i = 0; i<sw . t ; i++)
	{
		for( INT j=0; j<l; j++ )
		{
			if( thread_plcp[ i ][ j ] > PLCP[j] )	
			{
				PLCP[j] =  thread_plcp[ i ][ j ];
				P[j] = thread_p[i][j];
			}
		}
	}
	

	for(INT i =0; i<sw.t; i++)
	{
		free( thread_plcp[i] );
		free( thread_p[i] );
		free( column_lengths[i] );
		free( start_column[i] );
		free( order[i] );
		free( rank1[i] );
		free( rank2[i] );
		free( final_rank[i] );
		free( error_pos[i] );
	}
	free( thread_plcp );
	free( thread_p );
	free( column_lengths );
	free( start_column );
	free( order );
	free( rank1 );
	free( rank2 );
	free( final_rank );
	free( error_pos );
	delete( bucket );

return 1;
}
        


/*computing plcp for one set of error positions*/
INT compute_plcp( INT l, INT ** error_pos, INT * SA, INT * invSA, INT * LCP, INT * PLCP, INT * P, INT * A, INT ** thread_plcp, INT ** thread_p, TSwitch sw, INT ** column_lengths, INT ** start_column, INT ** order, INT ** rank1, INT ** rank2, INT ** final_rank, vector<vector<vector<INT>>> * bucket,  unordered_map<pair<INT, INT>, INT, pair_hash> * h_map   )
{

	INT thread = omp_get_thread_num();

	INT current_pos = 0 ; 
	INT plcp = 0;
	INT errorPosCur = - 1; 
	INT errorPosNxt = error_pos[thread][0];

	bool draw = false;
	bool first = true;

	for(INT b = 0; b<l; b++)
		rank1[thread][b] = 0;

	INT a = 0;

	for(a=0; a<=sw.k; a++)
	{
		INT column_length = (  errorPosNxt - errorPosCur ) - 1;

		column_lengths[thread][a] = column_length;
		start_column[thread][a] = current_pos;

		for(INT b = 0; b<l; b++)
		{
			order[thread][b] = -1;
			rank2[thread][b] = 0;
		}

		if( column_length > 0 )
		{			

			for(INT b = 0; b<l; b++)
			{	
				if( SA[b]+ current_pos  > l-1 )
					continue;
				else order[thread][invSA[SA[b]+current_pos]] = b;
				
			}
			
			INT rank = 1;

			if( first == true ) // first time need to add everything to add rank1 else add everything to rank2
			{
				for(INT b = 1; b<l; b++)
				{
					
					if( order[thread][b-1] != -1 )
					{
						if( LCP[b] >= column_length )
						{	
							rank1[thread][order[thread][b-1]] = rank;
							if( b > 1 )
								draw = true;
						}
						else rank1[thread][order[thread][b-1]] = rank++;
					}
				}

				if( order[thread][l-1] != -1)
				{
					if( LCP[l-1] >= column_length )
					{
						rank1[thread][order[thread][l-1]] = rank;
						draw = true;
					}
					else rank1[thread][order[thread][l-1]] = rank++;
				}

				first = false;

				if( draw == false )	
				{
					a++;
					break;
				}
			}
			else
			{
				
				for(INT b = 1; b<l; b++)
				{	
					if( order[thread][b-1] != -1 )
					{	
						if( LCP[b] >= column_length )
						{	
							rank2[thread][order[thread][b-1]] = rank;

							if( b > 1 )
								draw = true;
						}
						else rank2[thread][order[thread][b-1]] = rank++;
					}
				}

				if( order[thread][l-1] != -1)
				{
					if( LCP[l-1] >= column_length )
					{
						rank2[thread][order[thread][l-1]] = rank;
						draw = true;
					}
					else rank2[thread][order[thread][l-1]] = rank++;
				}

				if( draw == false )
				{
					a++;
					break;
				}
					
				bucket_sort_pair( rank1, rank2, l, order, final_rank, bucket, sw );
			}

			
		}

		errorPosCur = errorPosNxt;
		current_pos = errorPosCur + 1;
		if( a <sw.k-1 )
			errorPosNxt = error_pos[thread][a+1];
		else errorPosNxt = l;	
	}


	for(INT b = 0; b<l; b++)
	{	
		rank2[thread][b] = -1;
	}

	for(INT b = 0; b<l; b++)
	{	
		if( rank1[thread][b]-1 < 0 )
			continue ;
		else rank2[thread][rank1[thread][b]-1] = SA[b];
	}

	for(INT c = 1; c<l; c++ )
	{
		
		if( rank2[thread][c] == -1 || rank2[thread][c-1] == -1 )
			continue;

		INT total = 0;

		for(INT b = 0; b<a; b++) // k + 1 * n rmqs, break at a if no draws
		{	
			if( column_lengths[thread][b] == 0 )
				continue;
			
			INT rq = 0;


			if( rank2[thread][c-1]+start_column[thread][b] >= l || rank2[thread][c]+start_column[thread][b] >=l )
				rq = 0;
			else rq = LCP[ rmq(A, LCP, l, invSA[rank2[thread][c]+start_column[thread][b]], invSA[rank2[thread][c-1]+ start_column[thread][b]] ) ] ; 
			if(  rq >= column_lengths[thread][b] )
				total = total + column_lengths[thread][b];
			else if( rq < column_lengths[thread][b] )
			{
				total = total + rq ;
				break; // rq is not equal to whole column so error is occurring at end of column
			}
		}

		if( total + sw.k  >= sw.m )
		{
			pair<INT, INT> toadd( rank2[thread][c], rank2[thread][c-1] );
			unordered_map<pair<INT,INT>, INT, pair_hash>::const_iterator pos = h_map->find( toadd );
			if( pos == h_map->end()  )
			{
				pair< pair<INT, INT>, INT > toadd2( toadd, 1 );
				h_map->insert( toadd2 ); 
			}
		}

		if( min( total + sw.k, l - rank2[thread][c]  ) > thread_plcp[ thread  ][ rank2[thread][c] ] )
		{
			thread_plcp[ thread ][ rank2[thread][c] ] =  min( total + sw.k, l - rank2[thread][c] );	
			thread_p[ omp_get_thread_num() ][ rank2[thread][c] ] = rank2[thread][c-1];
		}

		if( min( total + sw.k, l - rank2[thread][c-1]  ) > thread_plcp[ thread ][ rank2[thread][c-1] ] )
		{

			thread_plcp[thread ][ rank2[thread][c-1] ] = min( total + sw.k, l - rank2[thread][c-1] );	
			thread_p[ thread ][ rank2[thread][c-1] ] = rank2[thread][c];
		}
	}

return 1;
}


INT bucket_sort_pair( INT ** rank1, INT ** rank2, INT l, INT ** order, INT ** final_rank, vector<vector<vector<INT>>> * bucket, TSwitch sw )
{

	INT thread = omp_get_thread_num();

	for(INT a = 0; a<l; a++)
	{
		final_rank[thread][a] = 0;
		bucket->at(thread).at(rank2[thread][a] ).push_back( a ) ;
	}
	

	INT a = 0;

	for(INT b = 0; b<l; b++)
	{
		for( INT c = 0; c<bucket->at(thread).at(b).size(); c++)
		{	
			order[thread][a] = bucket->at(thread).at(b).at(c);
			a++;
		}

		bucket->at(thread).at(b).clear();

	}

	for(INT b = 0; b<l; b++)
		bucket->at(thread).at( rank1[thread][order[thread][b]] ).push_back( order[thread][ b ] ) ;


	a = 0;

	for(INT b = 0; b<l; b++)
	{
		for( INT c = 0; c<bucket->at(thread).at(b).size(); c++)
		{
			order[thread][a] = bucket->at(thread).at(b).at(c);
			a++;
		}

		bucket->at(thread).at(b).clear();
	}

	INT rank = 0;

	for(INT b = 1; b<l; b++)
	{
		if( rank1[thread][order[thread][b]] == rank1[thread][order[thread][b-1]] && rank2[thread][order[thread][b]] == rank2[thread][order[thread][b-1]] )
		{
			final_rank[thread][order[thread][b]] = rank;
		}
		else 
		{
			rank = rank+1;
			final_rank[thread][order[thread][b]] = rank ;
		}
	}

	for(INT b = 0; b<l; b++)
	{	
		rank1[thread][b] = final_rank[thread][b];
	}


return 1;
}


INT rmqs( INT * LCP, INT first, INT second, INT l  )
{
	INT rmq = l;

	INT a = min( first, second);
	INT b = max( first, second );

	for(INT i =a+1; i<=b; i++)
	{
		if( LCP[ i ] < rmq )
			rmq = LCP[ i ];
	}

return rmq;
}

