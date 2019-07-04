# How-to:

To call the program you need python 2.7 and call opal.py with an input file *.opal 
For example: python opal.py test.opal

To run the software it is imperature that the spud library contained in the folder spud is compiled and the bindings for python a created. For this the following commands need to be run:
1) Open a terminal inside spud folder
2) Run: sudo ./configure
3) Run: sudo ./make 
4) Run: cd ./python
5) Run: sudo python setup.py build 
6) Run: sudo python setup.py install

This will create the necessary library for Opal to work properly.






