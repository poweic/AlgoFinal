#ifndef _CMD_PARSER_H
#define _CMD_PARSER_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

class CmdParser {

  private:

    class AutoType {
      public:

	AutoType(string str): _str(str) {}
	operator std::string () { return _str; }
	operator int    () { return str2int(_str); }
	operator size_t () { return str2int(_str); }
	operator float  () { return str2float(_str); }
	operator double () { return str2double(_str); }
	operator bool   () { 
	  assert(_str == "true" || _str == "false");
	  return _str == "true";
	}

      private:
	string _str;
    };
  
  public:
    CmdParser(int argc, char** argv): _argc(argc), _argv(argv), _usage(""), _options("") {
      _usage = "Usage: " + string(_argv[0]) + " ";
    }

    CmdParser& add(string option, bool mandatory = true) {
      static size_t counter = 0;
      string description = option;
      option = int2str(++counter);
      _appendUsage(!mandatory, description);

      Arg arg(option, description, !mandatory, "");
      _arguments[option] = arg;

      return *this;
    }

    CmdParser& add(string option, const char* description, bool mandatory = true) {
      bool optional = !mandatory;
      Arg arg(option, description, optional, "");
      _arguments[option] = arg;
      _options += arg.getDescription();
      _appendUsage(optional, option);

      return *this;
    }

    CmdParser& add(string option, const char* description, const char* defaultArg) {

      Arg arg(option, description, true, defaultArg);
      _arguments[option] = arg;
      _options += arg.getDescription();
      _appendUsage(true, option);

      return *this;
    }

    CmdParser& addGroup(string description) {
      if (description[description.size() - 1] == ':')
	description = description.substr(0, description.size() - 1);

      _options += "\n" + description + "\n";
      return *this;
    }

    AutoType operator[] (int argNumber) const {
      return find(int2str(argNumber));

    }

    AutoType operator[] (string option) const {
      return find(option);
    }

    AutoType find(string option) const {
      map<string, Arg>::const_iterator itr = _arguments.find(option);
      assert(itr != _arguments.end());
      return itr->second.parameter;
    }

    void showUsageAndExit() const {
      cout << endl;
      cout << _usage << endl;
      cout << _options << endl;
      exit(-1);
    }

    bool isOptionLegal() {
      if ( this->_lookingForHelp() )
	this->showUsageAndExit();

      return this->parse();
    }
    
    void showAll() {
      for (map<string, Arg>::iterator i=_arguments.begin(); i != _arguments.end(); ++i) {
	const string& opt = i->first;
	const string& parm = i->second.parameter;

	if (isNumber(opt))
	  printf("argument #%d:	%s\n", str2int(opt), parm.c_str());
	else
	  printf("option %s: %s\n", opt.c_str(), parm.c_str());
      }
    }

  private:

    struct Arg {
      Arg() {}

      Arg(string opt, string des, bool o, string darg) {
	option = opt;
	optional = o;
	default_arg = darg;
	description = replace_all(des, "\n", "\n" + getPadding());
	parameter = default_arg;
      }

      string _optionStr() const {
	string opt = "  " + option;
	opt.resize(PADDING_RIGHT, ' ');
	opt += "\t";
	return opt;
      }

      string getPadding() const {
	string padding;
	padding.resize(PADDING_RIGHT, ' ');
	padding += "\t";
	return padding;
      }

      string _defaultArgStr() const {
	if (optional && default_arg != "")
	  return getPadding() + "(default = " + default_arg + ")\n";

	return "";
      }

      string getDescription() const {
	return _optionStr() + description + "\n" + _defaultArgStr();
      }

      string getDefaultArg() const {
	return default_arg;
      }

      static const size_t PADDING_RIGHT = 24;

      string option;
      string description;
      bool optional;
      string default_arg;
      string parameter;
    };

    bool parse() {

      size_t counter = 0;

      for(int i=1; i<_argc; ++i) {
	string arg(_argv[i]);

	size_t pos = arg.find_first_of('=');
	if (pos == string::npos) {

	  if (has(arg)) {
	    if (i+1 >= _argc)
	      return miss(_arguments[arg]);
	    _arguments[arg].parameter = _argv[++i];
	  }
	  else {
	    if (!has(counter + 1))
	      return unknown(arg);

	    _arguments[int2str(counter + 1)].parameter = arg;
	    ++counter;
	  }
	}
	else {
	  string left = arg.substr(0, pos);
	  string right = arg.substr(pos + 1);

	  if (!has(left))
	    return unknown(left);

	  if (right.empty())
	    return miss(_arguments[left]);

	  _arguments[left].parameter = right;
	}
      }

      map<string, Arg>::iterator itr = _arguments.begin();
      for (; itr != _arguments.end(); ++itr) {
	const Arg& arg = itr->second;
	if (arg.parameter.empty() && !arg.optional)
	  return miss(arg);
      }

      return true;
    }

    bool miss(const Arg& arg) const {
      if (isNumber(arg.option))
	cerr << "Missing argument " + arg.description << endl;
      else
	cerr << "Missing argument after " + arg.option << endl;
      return false;
    }

    bool unknown(string opt) {
      cerr << "Unknown option " + opt << endl;
      return false;
    }

    bool has(string arg) {
      return _arguments.count(arg) > 0;
    }

    bool has(int n) {
      return _arguments.count(int2str(n)) > 0;
    }

    bool _lookingForHelp() const {
      string help("--help");
      for (int i=1; i<_argc; ++i) {
	if (help.compare(_argv[i]) == 0)
	  return true;
      }
      return false;
    }

    void _appendUsage(bool optional, string option) {
      if (optional)
	_usage += " [" + option + "]";
      else
	_usage += " <" + option + ">";
    }

    // ===== Static Utility Functions =====
    static bool isNumber(const std::string& s) {
      string::const_iterator it = s.begin();
      while (it != s.end() && std::isdigit(*it))
	++it;
      return !s.empty() && it == s.end();
    }

    static bool isFloat( string myString ) {
      std::istringstream iss(myString);
      float f;
      iss >> noskipws >> f; // noskipws considers leading whitespace invalid
      // Check the entire string was consumed and if either failbit or badbit is set
      return iss.eof() && !iss.fail(); 
    }

    static string int2str(int n) {
      stringstream ss;
      ss << n;
      return ss.str();
    }

    static string replace_all(const string& str, const string &token, const string &s) {
      string result(str);
      size_t pos = 0;
      while((pos = result.find(token, pos)) != string::npos) {
	result.replace(pos, token.size(), s);
	pos += s.size();
      }
      return result;
    }
    
    static int str2int(string str) {
      assert(isNumber(str));
      return ::atoi(str.c_str());
    }

    static float str2float(string str) {
      assert(isFloat(str));
      return ::atof(str.c_str());
    }

    static double str2double(string str) {
      assert(isFloat(str));
      return ::atof(str.c_str());
    }

    int _argc;
    char** _argv;

    string _usage;
    string _options;
    map<string, Arg> _arguments;
};

#endif // _CMD_PARSER_H
