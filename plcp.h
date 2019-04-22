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

#include <vector>
#include <math.h>
#include <unordered_set>

#define ALLOC_SIZE              1048576
 
#define MAX2(a,b) ((a) > (b)) ? (a) : (b)
#define MIN2(a,b) ((a) < (b)) ? (a) : (b)

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

struct TSwitch
{
    char *               input_filename;    
    char *               output_filename; 
    INT         	 k, t, m, r, l;
    double e;
};

struct Repeats
{
    INT pos1;    
    INT pos2;
    INT len;
    INT del;
};

struct Long
{
    INT pos1;    
    INT pos2;
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
	auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;  
    }
};



//INT rmqs( INT * LCP, INT first, INT second, INT l );

INT nchoosek( INT n, INT k );

double gettime( void );

int decode_switches ( int argc, char * argv [], struct TSwitch * sw  );

void usage( void );

INT LCParray( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );

INT compute_SA( unsigned char *text, INT n, INT * SA );

INT compute_invSA( unsigned char *text, INT n, INT * SA, INT * invSA );

INT compute_LCP( unsigned char *text, INT n, INT * SA, INT * invSA, INT * LCP );

INT populate_PLCP( unsigned char *text, INT l, INT * SA, INT * invSA, INT * LCP, INT * PLCP, INT * P );

INT LCParray( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );

INT repeats( INT * SA, INT * invSA, INT * LCP, INT * PLCP, INT * P, std::vector<Repeats> * sm_repeats, INT ell , TSwitch sw );

INT long_plcp( unsigned char * x, TSwitch  sw, INT * PLCP, INT * P, INT * SA, INT * LCP, std::unordered_set< std::pair<INT,INT> , pair_hash > * h_map  );

INT compute_long( INT * PLCP, INT * P, TSwitch sw, unsigned char * x, INT first, INT second )  ;

INT short_plcp( unsigned char * x, TSwitch sw, INT * PLCP, INT * P, INT * SA, INT * LCP, INT * invSA, INT * A, std::unordered_set< std::pair<INT,INT> , pair_hash > * h_map );

INT compute_plcp( INT l, INT ** error_pos, INT * SA, INT * invSA, INT * LCP, INT * PLCP, INT * P, INT * A, INT ** thread_plcp, INT ** thread_p, TSwitch sw, INT ** column_lengths, INT ** start_column, INT ** order, INT ** rank1, INT ** rank2, INT ** final_rank,  std::vector<std::vector<std::vector<INT>>> * bucket,  std::unordered_set< std::pair<INT,INT> , pair_hash > * h_map   );

INT bucket_sort_pair(INT ** rank1, INT ** rank2, INT l, INT ** order, INT ** final_rank,  std::vector<std::vector<std::vector<INT>>> * bucket, TSwitch sw );



//INT prefix_len( unsigned char * x, struct TSwitch sw, INT first, INT second );

static __inline void  binary_search ( INT * left, INT * right, INT * key, INT value )
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

#ifndef RMQ_H
#define RMQ_H

void rmq_preprocess(INT *m, INT *v, INT n);
//INT rmq(INT *m, INT *v, INT n, INT i, INT j);

static __inline INT flog2(INT v) {
  return (INT) floor (log2((double)(v)));
}

static __inline INT rmq(INT *m, INT *v, INT n, INT i, INT j) {
  INT lgn = flog2(n); 
  if (i > j) {INT tmp = j; j = i; i = tmp;}
  i++; //for LCP
  if (i == j) return i;
  INT k = flog2(j-i+1); 
  INT a = m[i * lgn + k]; 
  INT shift = ((( INT ) 1) << k);
  INT b = m[(j - shift + 1) * lgn + k];  
  return v[a]>v[b]?b:a;
}

#endif

