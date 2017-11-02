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
    unsigned int         k, m, T;
};

double gettime( void );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw  );
void usage ( void );
unsigned int LCParray( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );

unsigned int compute_SA( unsigned char *text, INT n, INT * SA );

unsigned int compute_invSA( unsigned char *text, INT n, INT * SA, INT * invSA );

unsigned int compute_LCP( unsigned char *text, INT n, INT * SA, INT * invSA, INT * LCP );

unsigned int populate_PLCP( unsigned char *text, int l, INT * SA, INT * invSA, INT * LCP, unsigned int * PLCP, unsigned int * P );

unsigned int LCParray( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );

unsigned int k_mappability( unsigned char * x, struct TSwitch  sw, unsigned int * PLCP, unsigned int * P, INT * SA, INT * LCP );

unsigned int short_plcp( int i, char * alphabet, unsigned char * x, struct TSwitch  sw, unsigned int * PLCP, unsigned int * P, INT * SA, INT * LCP, INT * invSA, INT * A);

unsigned int extension( char * alphabet, unsigned int * character, unsigned int * error_pos,  int i , INT * SA, INT * invSA, INT * LCP, unsigned int * PLCP, unsigned int * P, struct TSwitch  sw, unsigned char * x,  int error, unsigned int * prev_errors, INT * A);

unsigned int positions( char * alphabet, unsigned int * character, unsigned int * error_pos, int i, int error,  INT *  SA, INT * invSA, INT * LCP, unsigned int * PLCP, unsigned int * P, struct TSwitch sw, unsigned char * x, int start, int  end, int len_ext, unsigned int * prev_errors,  INT * A );

#ifndef RMQ_H
#define RMQ_H

void rmq_preprocess(INT *m, INT *v, INT n);
INT rmq(INT *m, INT *v, INT n, INT i, INT j);

static __inline INT flog2(INT v) {
  return (INT) floor (log2((double)(v)));
}

#endif

