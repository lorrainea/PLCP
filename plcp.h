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

#define ALLOC_SIZE              1048576

#define DNA                     "ACGT"                        
#define PROT                    "ARNDCQEGHILKMFPSTWYVOUBZJX*"     
#define MAX2(a,b) ((a) > (b)) ? (a) : (b)
#define MIN2(a,b) ((a) < (b)) ? (a) : (b)

struct TSwitch
{
    char *               input_filename;         // the input file name
    char *               output_filename;        // the output file name
    char *	         alphabet;
    unsigned int         k, m;
};

double gettime( void );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw  );
void usage ( void );

