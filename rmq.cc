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

#include <stdlib.h>
#include <assert.h>
#include "plcp.h"
#include <stdio.h>

void rmq_preprocess( INT * m, INT * v, INT n ) 
{
	INT i, j;
	INT lgn = flog2(n); 
	for (i = 0; i < n; i++)
		m[i*lgn] = i;
	for (j = 1; ((( INT ) 1) << j) <= n; j++)
		for (i = 0; i + (1 << j) - 1 < n; i++)
			if (v[m[i*lgn + j - 1]] < v[m[(i + (1 << (j - 1)))*lgn + j - 1]])
				m[i*lgn + j] = m[i*lgn + j - 1];
			else
				m[i*lgn + j] = m[(i + (1 << (j - 1)))*lgn + j - 1];
}
