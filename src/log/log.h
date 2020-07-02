//------------------------------------------------------------------------------
//  log.h
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara, F. Costa, P. Pessoa]
//  [ncarrara@albany.edu,felipecosta.physics@gmail.com,
//    pedroh.pessoa100@gmail.com]
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
#pragma once

//------------------------------------------------------------------------------
//	The basic design of this logging system was taken from
//	Cherno's Hazel engine: https://github.com/TheCherno/Hazel/blob/master/Hazel/src/Hazel/Core/Log.h
//------------------------------------------------------------------------------

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <complex>
#include <lapacke.h>
#include <cblas.h>
#include <stdint.h>
#include <complex>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
#include <sys/stat.h>
#include <sys/types.h>

namespace ET
{

  class Log
	{
	public:
		Log();
		~Log();
		//--------------------------------------------------------------------------
    //  initializer for the logging system for an object
		//--------------------------------------------------------------------------
    void init(std::string name, const std::string& outputFile=".log/default.txt",
							uint32_t debug=2);

		//--------------------------------------------------------------------------
    //  getter for the pointer for the logger object
		//--------------------------------------------------------------------------
    std::shared_ptr<spdlog::logger>& getLogger() { return _logger; }
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Generators for logging messages
		//--------------------------------------------------------------------------
		void TRACE(std::string data);
		void INFO(std::string data);
		void WARN(std::string data);
		void ERROR(std::string data);
		void CRITICAL(std::string data);
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Get the last 'n' lines of the outputFile
		//--------------------------------------------------------------------------
		std::string getOutput(uint64_t lines);
		std::string getOutput();
		//--------------------------------------------------------------------------

  private:
		std::string _name;
		std::string _outputFile;
    //  shared pointer to our logger object
    std::shared_ptr<spdlog::logger> _logger;
	};

}
