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
#include <math.h>
#define ALLOC_SIZE              1048576

#define DNA                     "ACGT"                        
#define PROT                    "ARNDCQEGHILKMFPSTWYVOUBZJX*"     
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
    char *               input_filename;         // the input file name
    char *               output_filename;        // the output file name
//    char *	         alphabet;
    INT         k, m, T;
};

double gettime( void );

int decode_switches ( int argc, char * argv [], struct TSwitch * sw  );

void usage ( void );

INT LCParray( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );

INT compute_SA( unsigned char *text, INT n, INT * SA );

INT compute_invSA( unsigned char *text, INT n, INT * SA, INT * invSA );

INT compute_LCP( unsigned char *text, INT n, INT * SA, INT * invSA, INT * LCP );

INT populate_PLCP( unsigned char *text, INT l, INT * SA, INT * invSA, INT * LCP, INT * PLCP, INT * P );

INT LCParray( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );

INT k_mappability( unsigned char * x, struct TSwitch  sw, INT * PLCP, INT * P, INT * SA, INT * LCP );

INT short_plcp( INT i, char * alphabet, unsigned char * x, struct TSwitch  sw, INT * PLCP, INT * P, INT * SA, INT * LCP, INT * invSA, INT * A);

INT extension( char * alphabet, INT * character, INT * error_pos,  INT i , INT * SA, INT * invSA, INT * LCP, INT * PLCP, INT * P, struct TSwitch  sw, unsigned char * x,  INT error, INT * prev_errors, INT * A);

INT positions( char * alphabet, INT * character, INT * error_pos, INT i, INT error,  INT *  SA, INT * invSA, INT * LCP, INT * PLCP, INT * P, struct TSwitch sw, unsigned char * x, INT start, INT end, INT len_ext, INT * prev_errors,  INT * A );

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

