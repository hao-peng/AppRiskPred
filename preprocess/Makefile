all: ../bin/split ../bin/changeFormat ../bin/combine

../bin/combine: combine.cpp
	g++ -o ../bin/combine combine.cpp

../bin/split: split.cpp
	g++ -o ../bin/split split.cpp

../bin/changeFormat: changeFormat.cpp
	g++ -o ../bin/changeFormat changeFormat.cpp

clean:
	rm -f ../bin/split
	rm -f ../bin/changeFormat
	rm -f ../bin/combine


