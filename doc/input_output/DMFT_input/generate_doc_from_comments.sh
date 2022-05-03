#!/bin/bash

# Changes the input from:
#


sed 's/\(.*\) :\s*\(.*\)/\n.. admonition:: \1: \n \n            \*\*type=\*\* \2\n/g' python_comments.txt > matches_comments.txt

#add blank line after type
#sed -i 's/\(\*type=.*\),/\1 \n\n           /g'  matches_comments.txt
#sed -i 's/\(\*type=.*\),/\1 \n\n           /g'  matches_comments.txt

# make 'optional' and 'default' bold
sed -i 's/,.*\(\optional\)/;  \*\*\1\*\*/g'  matches_comments.txt
sed -i 's/,.*\(\default=\)/;  \*\*\1\*\* /g'  matches_comments.txt

# grep all admonitions and store them in a file 

cat > input.rst << EOF
DMFT input
------------------------

In this section all the 

Input/Output
===================
.. toctree::
    :maxdepth: 1

    general
    solver
    dft
    advanced

Below an exhaustive list containing all the parameters marked by section.

**[general]**

EOF
# script to generate subsection divided by the '[  SECTION  ]' pattern
#
# awk '/\[  SECTION  \]/{flag=1; c=0} flag; /\[ /&& ++c==2{flag=0}'
# matches from the pattern '[  SECTION  ]'(note, two spaces are important) up to either the next occurence of
# '[  ', which is the next group, or the end of the file 
#
#
# the second part takes the divided section and extracts a list out of it:
# grep '::', to match all annotations
# grep -o '::\(.*\):', to only match the part after the '::'
# sed 's/:: //g' to remove the ::
# tr ':\n' '; ' to remove the newline and add a semicolon

awk '/\[  general  \]/{flag=1; c=0} flag; /\[ /&& ++c==2{flag=0}' matches_comments.txt | head -n -2 | tail -n +3 > general.tmp
grep '::' general.tmp | grep -o '::\(.*\):' | sed 's/:: //g' | tr ':\n' '; ' > general_list.tmp
awk '/\[  solver  \]/{flag=1; c=0} flag; /\[ /&& ++c==2{flag=0}' matches_comments.txt | head -n -2 | tail -n +3 > solver.tmp
grep '::' general.tmp | grep -o '::\(.*\):' | sed 's/:: //g' | tr ':\n' '; ' > solver_list.tmp
awk '/\[  dft  \]/{flag=1; c=0} flag; /\[ /&& ++c==2{flag=0}' matches_comments.txt | head -n -2 | tail -n +3 > dft.tmp
grep '::' general.tmp | grep -o '::\(.*\):' | sed 's/:: //g' | tr ':\n' '; ' > dft_list.tmp
awk '/\[  advanced  \]/{flag=1; c=0} flag; /\[ /&& ++c==2{flag=0}' matches_comments.txt | head -n -2 | tail -n +3 > advanced.tmp
grep '::' general.tmp | grep -o '::\(.*\):' | sed 's/:: //g' | tr ':\n' '; ' > advanced_list.tmp


cat general_list.tmp >> input.rst

echo -e "\n"  >> input.rst
echo  '**[solver]**' >> input.rst
echo -e "\n" >> input.rst
cat solver_list.tmp >> input.rst

echo -e "\n"  >> input.rst
printf '**[dft]**' >> input.rst
echo -e "\n"  >> input.rst
cat dft_list.tmp >> input.rst

echo -e "\n"  >> input.rst
printf '**[advanced]**' >> input.rst
echo -e "\n"  >> input.rst
cat advanced_list.tmp >> input.rst

#cat matches_comments.txt >> input.rst

# now create the subpages, one for every group in the ini file
# there are two awk scripts, the second one is the simplest
#

##############
cat > general.rst << EOF
[general]: General parameters
------

Includes the majority of the parameters

List of possible entries:

EOF
cat general_list.tmp >> general.rst
echo -e "\n"  >> general.rst
cat general.tmp >> general.rst
##############

##############
cat > solver.rst << EOF
[solver]: solver specific parameters
------

Here are the parameters that are uniquely dependent on the solver chosen. Below a list of the supported solvers:

===================
.. toctree::
    :maxdepth: 1



List of possible entries:

EOF
cat solver_list.tmp >> solver.rst
echo -e "\n"  >> solver.rst
cat solver.tmp >> solver.rst
##############



###############
cat > dft.rst << EOF

[dft]: DFT related inputs
------

List of parameters that relate to the DFT calculation, useful mostly when doing CSC.

List of possible entries:


EOF

cat dft_list.tmp >> dft.rst
echo -e "\n"  >> dft.rst
cat dft.tmp >> dft.rst
##############

##############
cat > advanced.rst << EOF
[advanced]: Advanced inputs
------

Advanced parameters, do not modify default value unless you know what you are doing

List of possible entries:

EOF
cat advanced_list.tmp >> advanced.rst
echo -e "\n"  >> advanced.rst
cat advanced.tmp >> advanced.rst
##############

rm ./*.tmp
