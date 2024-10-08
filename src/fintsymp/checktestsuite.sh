#! /bin/bash
#pour tester uniquement avec gfortran
alltest="0"
echo $FC | grep gfortran && alltest="1"
echo $alltest

#execution des tests

cd intplabasic && make clean && make && cd ..

cd intplastat && make clean && make && make check  && cd .. 

cd intpartbasic && make clean && make && make check && cd .. 

cd intpartstat && make clean && make  && cd .. #&& make check 
if test "x$alltest" == "x1" ; then
cd intpartstat && make check && cd ..
fi

cd intcircumbinbasic && make clean && make && make check && cd .. 

cd intcircumbinmaree && make clean && make && cd ..

cd heliodebasic && make clean && make && cd ..
