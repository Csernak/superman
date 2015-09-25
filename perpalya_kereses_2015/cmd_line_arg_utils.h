//============================================================================
// Name        : cycle.cpp
// Author      : Gergely Gyebrószki
// Version     :
// Copyright   : 2014 (C) gyebro.com
// Description : Simple Command Line Argument Parser
//============================================================================

#ifndef CMD_LINE_ARG_UTILS_H
#define CMD_LINE_ARG_UTILS_H

#include <vector>
#include <string>
#include <iostream>
#include "stdlib.h"

using namespace std;

class CommandLineArgumentParser {
private:
	int argcount;
	vector<string> arguments;
public:
	CommandLineArgumentParser(int argc, wchar_t** argv) {
		argcount = argc;
		arguments.clear();
		arguments.resize(0);
		for (int i=0; i<argc; i++) {
			wstring ws(argv[i]);
			string s(ws.begin(),ws.end());
			arguments.push_back(s);
		}
	}
	CommandLineArgumentParser(int argc, char** argv) {
		argcount = argc;
		arguments.clear();
		arguments.resize(0);
		for (int i=0; i<argc; i++) {
			string s(argv[i]);
			arguments.push_back(s);
		}
	}
	bool GetString(const string keyword, string& s) {
		bool ok = false;
		for (int i=1; i<argcount-1; i++) {
			if (arguments[i] == keyword) {
				s = arguments[i+1];
				ok = true;
			}
		}
		return ok;
	}
};

#endif
