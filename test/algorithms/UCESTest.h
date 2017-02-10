//
// SFFSTest.h -- definition of the namespace "SFFSTest".
//
//    This file is part of the featsel program
//    Copyright (C) 2010  Marcelo S. Reis
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef UCESTEST_H_
#define UCESTEST_H_

#include "../../src/algorithms/UCES.h"
#include "../../src/functions/SubsetSum.h"
#include "../../src/functions/HammingDistance.h"
#include "../../src/functions/Explicit.h"

namespace UCESTest
{

	bool check_number_of_minima();

    bool check_results();

    bool check_results2();

    bool check_results3();
}

#endif /* UCESTEST_H_ */
