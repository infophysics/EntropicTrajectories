//------------------------------------------------------------------------------
//  log.cpp
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

	void Log::init(std::string name, const std::string& outputFile)
	{
		_name = name;
		std::vector<spdlog::sink_ptr> logSinks;
		logSinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
		logSinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(outputFile, true));

		logSinks[0]->set_pattern("%^[%T] %n: %v%$");
		logSinks[1]->set_pattern("[%T] [%l] %n: %v");

		_logger = std::make_shared<spdlog::logger>(_name, begin(logSinks), end(logSinks));
		spdlog::register_logger(_logger);
		_logger->set_level(spdlog::level::trace);
		_logger->flush_on(spdlog::level::trace);
	}

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
}
