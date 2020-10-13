//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>
#include <algorithm>
#include <sstream>
#include <execinfo.h>

#include "logger.h"

template<class T>
void vecRemove(std::vector<T>& v, T val)
{
	v.erase(std::remove(v.begin(), v.end(), val), v.end()); 
}

struct pairhash 
{
public:
	template <typename T, typename U>
	std::size_t operator()(const std::pair<T, U> &x) const
	{
		return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
	}
};


template<typename T>
T quantile(const std::vector<T>& vec, int percent)
{
	if (vec.empty()) return 0;
	//NOTE: there's a bug in libstdc++ nth_element, 
	//that sometimes leads to a segfault. This is why
	//we have this inefficient impleemntation here
	//std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, 
	//				 vec.end());
	auto sortedVec = vec;
	std::sort(sortedVec.begin(), sortedVec.end());
	size_t targetId = std::min(vec.size() * (size_t)percent / 100, 
							   vec.size() - 1);
	return sortedVec[targetId];
}

template<typename T>
T median(const std::vector<T>& vec)
{
	return quantile(vec, 50);
}

inline std::vector<std::string> 
splitString(const std::string &s, char delim) 
{
	std::vector<std::string> elems;
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) elems.push_back(item);
	return elems;
}

inline bool fileExists(const std::string& path)
{
	std::ifstream fin(path);
	return fin.good();
}

inline void segfaultHandler(int signal __attribute__((unused)))
{
	void *stackArray[20];
	size_t size = backtrace(stackArray, 10);
	Logger::get().error() << "Segmentation fault! Backtrace:";
	char** backtrace = backtrace_symbols(stackArray, size);
	for (size_t i = 0; i < size; ++i)
	{
		Logger::get().error() << "\t" << backtrace[i];
	}
	abort();
}

inline void exceptionHandler()
{
	static bool triedThrow = false;
	try
	{
        if (!triedThrow)
		{
			triedThrow = true;
			throw;
		}
    }
    catch (const std::exception &e) 
	{
        Logger::get().error() << "Caught unhandled exception: " << e.what();
    }
	catch (...) {}

	void *stackArray[20];
	size_t size = backtrace(stackArray, 10);
	char** backtrace = backtrace_symbols(stackArray, size);
	for (size_t i = 0; i < size; ++i)
	{
		Logger::get().error() << "\t" << backtrace[i];
	}
	abort();
}

