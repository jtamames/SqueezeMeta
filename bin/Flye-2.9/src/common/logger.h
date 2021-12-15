//(c) 2013-2016 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <iostream>
#include <fstream>
#include <ctime>
#include <stdexcept>

class Logger
{
public:
	static Logger& get()
	{
		static Logger instance;
		return instance;
	}

	void setOutputFile(const std::string& filename)
	{
		_logFile.open(filename, std::ofstream::out | std::ofstream::app);
		_logFileSet = true;
		if (!_logFile.is_open()) 
		{
			throw std::runtime_error("Can't open log file");
		}
	}

	void setDebugging(bool debug) {_debug = debug;}

	class StreamWriter
	{
	public:
		StreamWriter(const std::string& level, 
					 std::ostream* consoleStream = nullptr,
					 std::ostream* fileStream = nullptr):
			_fileStream(fileStream), _consoleStream(consoleStream)
		{
			if (_fileStream) *_fileStream << timestamp() << " " << level << " ";
			if (_consoleStream) *_consoleStream << timestamp() 
												<< " " << level << " ";
		}
		~StreamWriter()
		{
			if (_fileStream) *_fileStream << std::endl;
			if (_consoleStream) *_consoleStream << std::endl;
		}
		
		template <class T>
		Logger::StreamWriter& operator<< (const T& val)
		{
			if (_fileStream) *_fileStream << val;
			if (_consoleStream) *_consoleStream << val;
			return *this;
		}

	private:
		std::ostream* _fileStream;
		std::ostream* _consoleStream;
	};

	StreamWriter info()
	{
		std::ostream* logPtr = _logFileSet ? &_logFile : nullptr;
		return StreamWriter("INFO:", &std::cerr, logPtr);
	}

	StreamWriter warning()
	{
		std::ostream* logPtr = _logFileSet ? &_logFile : nullptr;
		return StreamWriter("WARNING:", &std::cerr, logPtr);
	}

	StreamWriter error()
	{
		std::ostream* logPtr = _logFileSet ? &_logFile : nullptr;
		return StreamWriter("ERROR:", &std::cerr, logPtr);
	}

	StreamWriter debug()
	{
		std::ostream* logPtr = _logFileSet ? &_logFile : nullptr;
		std::ostream* consolePtr = _debug ? &std::cerr : nullptr;
		return StreamWriter("DEBUG:", consolePtr, logPtr);
	}

private:
	static std::string timestamp(const char* format = "[%Y-%m-%d %H:%M:%S]")
	{
		std::time_t t = std::time(0);
		char cstr[256];
		std::strftime(cstr, sizeof(cstr), format, std::localtime(&t));
		return cstr;
	}

	Logger():
		_debug(false), _logFileSet(false)
	{}
	~Logger()
	{
		if (_logFileSet)
		{
			_logFile << "-----------End assembly log------------\n";
		}
	}

	bool _debug;
	bool _logFileSet;
	std::ofstream _logFile;
};
