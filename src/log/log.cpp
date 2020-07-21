//------------------------------------------------------------------------------
//  log.cpp
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara]
//  [ncarrara@albany.edu]

//
//  Permission to use, copy, modify, and/or distribute this software for any
//  purpose with or without fee is hereby granted.
//
//  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
//  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
//  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
//  SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
//  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
//  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
//  IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//	The basic design of this logging system was taken from
//	Cherno's Hazel engine: https://github.com/TheCherno/Hazel/blob/master/Hazel/src/Hazel/Core/Log.cpp
//------------------------------------------------------------------------------

#include "log.h"

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

namespace ET
{

	Log::Log()
	{
	}
	Log::~Log()
	{
	}
	//----------------------------------------------------------------------------
	//  Initializer for the logger.
	//	the precompiler directive argument determines the verbosity of the logger.
	//	If one wants the full functionality then one must define 'LOG_LEVEL_DEBUG'
	//	during compilation, otherwise only an output file will be generated.
	//----------------------------------------------------------------------------
	void Log::init(std::string name, const std::string& outputFile, uint32_t debug)
	{
		_name = name;
		_outputFile = outputFile;
		std::vector<spdlog::sink_ptr> logSinks;
		//	Generate the file logger
		//	The mkdir command only works on Linux.  Would need to add support for
		//	Windows and MAC.
		if (mkdir(".logs",0777) == -1) {
			struct stat info;
			if (stat(".logs", &info) != 0) {
				std::cout << "ERROR: directory logs could not be created!" << std::endl;
			}
		}
		//	check that the logger doesn't exist already
		auto l = spdlog::get(name);
		//	if so, try and add numbers to the name
		std::string new_name = name;
		int i = 0;
		while (l != NULL) {
			new_name = name + "_" + std::to_string(i);
			i += 1;
			l = spdlog::get(new_name);
		}
		_name = new_name;
		logSinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>
			                    (outputFile, true));
#ifdef LOG_LEVEL_DEBUG
		logSinks[0]->set_pattern("[%T] [%l] %n: %v");
		logSinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
		logSinks[1]->set_pattern("%^[%T] %n: %v%$");
#else
		logSinks[0]->set_pattern("[%T] [%l] %n: %v");
#endif

		_logger = std::make_shared<spdlog::logger>(_name, begin(logSinks),
		                                           end(logSinks));
		spdlog::register_logger(_logger);
		_logger->set_level(spdlog::level::trace);
		_logger->flush_on(spdlog::level::trace);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Generators for logging messages
	//----------------------------------------------------------------------------
	void Log::TRACE(std::string data)
	{
		_logger->trace(data);
	}
	void Log::INFO(std::string data)
	{
		_logger->info(data);
	}
	void Log::WARN(std::string data)
	{
		_logger->warn(data);
	}
	void Log::ERROR(std::string data)
	{
		_logger->error(data);
	}
	void Log::CRITICAL(std::string data)
	{
		_logger->critical(data);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Get the last 'n' lines of the outputFile
	//----------------------------------------------------------------------------
	std::string Log::getOutput(uint64_t lines)
	{
    std::ifstream file(_outputFile);
		std::vector<std::string> l;
		std::string str;
		if (!file.good()) {
			return "ERROR: Could not open log file";
		}
		while (std::getline(file,str)) {
			l.push_back(str);
		}
		std::string result;
		if (lines > l.size()) {
			lines = l.size();
		}
		for (uint64_t i = l.size()-lines; i < l.size(); i++) {
			result += l[i] + "\n";
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Get entire output
	//----------------------------------------------------------------------------
	std::string Log::getOutput()
	{
		std::ifstream file(_outputFile);
		std::string result;
		std::string str;
		if (!file.good()) {
			return "ERROR: Could not open log file";
		}
		while (std::getline(file,str)) {
			result += str + "\n";
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//--------------------------------------------------------------------------
  //  Error messages
  //  Here we have a set of functions for generating various
  //  error messages
  //--------------------------------------------------------------------------
  std::string
	MATRIX_INCONSISTENT_ARRAY(std::vector<std::pair<size_t,size_t>>& rows)
  {
    std::string error = "ERROR: Attempted to construct matrix with";
    error += " inconsistent numbers of columns.";
    for (size_t i = 0; i < rows.size(); i++) {
      error += "  There are " + std::to_string(std::get<1>(rows[i]))
            +  " rows of size " + std::to_string(std::get<0>(rows[i])) + ".";
    }
    return error;
  }
	//--------------------------------------------------------------------------
  std::string
	MATRIX_OUT_OF_BOUNDS(bool axis, const size_t& bound, const size_t& attempt,
                       const std::string& name)
  {
    std::string error = "ERROR: Attempted to access index "
                      + std::to_string(attempt);
    if (name != " ") {
      error += " of matrix '" + name + "'";
    }
    error += "; out of bounds for axis " + std::to_string(int(axis))
           + " with size " + std::to_string(bound);
    return error;
  }
	//--------------------------------------------------------------------------
  std::string
	MATRIX_ADD_INCOMPATIBLE_ROWS(const size_t& m1, const size_t& m2,
                               const std::string& name1,
                               const std::string& name2)
  {
    std::string error = "ERROR: Attempted to add incompatible matrices";
    if (name1 != " " && name2 != " ") {
      error += "'" + name1 + "' and '" + name2 + "'";
    }
    error += "; row sizes are m1 = " + std::to_string(m1)
           + " and m2 = " + std::to_string(m2);
    return error;
  }
	//--------------------------------------------------------------------------
  std::string
	MATRIX_ADD_INCOMPATIBLE_COLS(const size_t& n1, const size_t& n2,
                             	 const std::string& name1,
                               const std::string& name2)
  {
    std::string error = "ERROR: Attempted to add incompatible matrices";
    if (name1 != " " && name2 != " ") {
      error += "'" + name1 + "' and '" + name2 + "'";
    }
    error += "; column sizes are n1 = " + std::to_string(n1)
           + " and n2 = " + std::to_string(n2);
    return error;
  }
	//--------------------------------------------------------------------------
  std::string
	MATRIX_SUB_INCOMPATIBLE_ROWS(const size_t& m1, const size_t& m2,
                        			 const std::string& name1,
                               const std::string& name2)
  {
    std::string error = "ERROR: Attempted to subtract incompatible matrices";
    if (name1 != " " && name2 != " ") {
      error += "'" + name1 + "' and '" + name2 + "'";
    }
    error += "; row sizes are m1 = " + std::to_string(m1)
           + " and m2 = " + std::to_string(m2);
    return error;
  }
	//--------------------------------------------------------------------------
  std::string
	MATRIX_SUB_INCOMPATIBLE_COLS(const size_t& n1, const size_t& n2,
                               const std::string& name1,
                               const std::string& name2)
  {
    std::string error = "ERROR: Attempted to subtract incompatible matrices";
    if (name1 != " " && name2 != " ") {
      error += "'" + name1 + "' and '" + name2 + "'";
    }
    error += "; column sizes are n1 = " + std::to_string(n1)
           + " and n2 = " + std::to_string(n2);
    return error;
  }
	//--------------------------------------------------------------------------
  std::string
	MATRIX_MUL_INCOMPATIBLE(const size_t& n1, const size_t& m2)
  {
    return " ";
  }
	//--------------------------------------------------------------------------
  std::string
	MATRIX_ZERO_DIV(const size_t& m, const size_t& n)
  {
    return " ";
  }
  //--------------------------------------------------------------------------
}
