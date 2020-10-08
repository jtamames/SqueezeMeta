//(c) 2020 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <string>

int assemble_main(int argc, char** argv);
int repeat_main(int argc, char** argv);
int contigger_main(int argc, char** argv);
int polisher_main(int argc, char** argv);

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::cerr << "Usage: flye-modules [assemble | repeat | contigger | polisher] ..." 
				  << std::endl;
		return 1;
	}
	std::string module = argv[1];
	if (module == "assemble")
	{
		return assemble_main(argc - 1, argv + 1);
	}
	else if (module == "repeat")
	{
		return repeat_main(argc - 1, argv + 1);
	}
	else if (module == "contigger")
	{
		return contigger_main(argc - 1, argv + 1);
	}
	else if (module == "polisher")
	{
		return polisher_main(argc - 1, argv + 1);
	}
	else
	{
		std::cerr << "Usage: flye-modules [assemble | repeat | contigger | polisher] ..." 
				  << std::endl;
		return 1;
	}

	return 0;
}
