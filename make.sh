if [[ -f coefficient.a ]]
then
    echo "build coefficient.a"
    gfortran -c reconstruction.f90
    ar rcs coefficient.a reconstruction.o
fi
c++ -c *cpp
c++ -o coeffgen *.o -lstdc++