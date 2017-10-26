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

#include <vector>

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

unsigned int LCParray( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );

unsigned int k_mappability( int i, unsigned char * x, struct TSwitch  sw, unsigned int * PLCP, unsigned int * P, INT * SA, INT * LCP );

unsigned int short_plcp( int i, char * alphabet, unsigned char * x, struct TSwitch  sw, unsigned int * PLCP, unsigned int * P, INT * SA, INT * LCP, INT * invSA);

unsigned int extension( char * alphabet, unsigned int * character, unsigned int * error_pos,  int i , INT * SA,unsigned int * PLCP, unsigned int * P, struct TSwitch  sw, unsigned char * x, int error);

unsigned int positions( char * alphabet, unsigned int * character, unsigned int * error_pos, int i, int error,  INT *  SA, unsigned int * PLCP, unsigned int * P, struct TSwitch sw, unsigned char * x, int start, int  end, int len_ext);
