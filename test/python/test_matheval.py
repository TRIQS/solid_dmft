# Copyright (c) 2018-2022 Simons Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You may obtain a copy of the License at
#     https:#www.gnu.org/licenses/gpl-3.0.txt
#
# Authors: Alexander Hampel

from solid_dmft.dmft_tools.matheval import MathExpr
import unittest


class test_mathexpr(unittest.TestCase):
    def test_simple(self):
        expr = MathExpr("1 + 1")
        result = expr()
        self.assertEqual(result, 2)

    def test_variables(self):
        expr = MathExpr("34788 * it + 928374 * rank")
        result = expr(it=5, rank=9)
        self.assertEqual(result, 34788 * 5 + 928374 * 9)

    def test_breakout(self):
        with self.assertRaises(ValueError):
            expr = MathExpr("(1).__class__")


if __name__ == "__main__":
    unittest.main()
