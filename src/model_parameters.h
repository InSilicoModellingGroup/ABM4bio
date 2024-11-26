// =============================================================================
//
//   Copyright (C) 2020-2024 Vasileios Vavourakis (vasvav@gmail.com)
//   All Rights Reserved.
//
//   Licensed under the GNU General Public License v3.0 (the "License").
//   See the LICENSE file provided in this project details the License.
//   You cannot use this file except in compliance with the License.
//
// =============================================================================

// =============================================================================
#ifndef _MODEL_PARAMETERS_H_
#define _MODEL_PARAMETERS_H_
// =============================================================================
#include <cstddef>
#include <string>
#include <vector>
#include <map>
#include <typeinfo>
#include <sstream>
#include "biodynamo.h"
// =============================================================================
// Helper functions for printing scalar, vector and vector<vector> types.
// Called from Parameters::Parameter<T>::print(...).
template<typename P>
void print_helper(std::ostream& os, const P* param);
template<typename P>
void print_helper(std::ostream& os, const std::vector<P>* param);
template<typename P>
void print_helper(std::ostream& os, const std::vector<std::vector<P>>* param);
// =============================================================================
// This class provides the ability to map between arbitrary, user-defined
// strings and several data types. This can be used to provide arbitrary
// user-specified options. Author: Benjamin S. Kirk (2004)
class Parameters {
//
public:
    // constructors & destructor
    Parameters() {}
    Parameters(const Parameters& p) { *this = p; }
    virtual ~Parameters() { clear(); }
    // Assignment operator.  Removes all parameters in \p this
    // and inserts copies of all parameters from \p source
    virtual Parameters& operator=(const Parameters & source);
    // Addition/Assignment operator.  Inserts copies of all parameters
    // from \p source.  Any parameters of the same name already in \p
    // this are replaced.
    virtual Parameters& operator+=(const Parameters & source);
    // \returns \p true if a parameter of type \p T
    // with a specified name exists, \p false otherwise.
    template <typename T>
    bool have_parameter(const std::string& ) const;
    // \returns A constant reference to the specified parameter
    // value. Requires, of course, that the parameter exists.
    template <typename T>
    const T& get(const std::string& ) const;
    // Inserts a new Parameter into the object but does not return
    // a writable reference.  The value of the newly inserted
    // parameter may not be valid.
    template <typename T>
    void insert(const std::string& );
    // \returns A writable reference to the specified parameter.
    // This method will create the parameter if it does not exist,
    // so it can be used to define parameters which will later be
    // accessed with the \p get() member.
    template <typename T>
    T& set(const std::string& );
    // Overridable function to set any extended attributes for
    // classes inheriting from this class.
    virtual void set_attributes(const std::string&, bool /*inserted_only*/) {}
    // Removes the specified parameter from the list, if it exists.
    void remove(const std::string& );
    // \returns The total number of parameters.
    std::size_t n_parameters() const { return _values.size(); }
    // Clears internal data structures & frees any allocated memory.
    virtual void clear();
    // Prints the contents, by default to std::cout.
    void print(std::ostream& os = std::cout) const;
    //*******************************************************
    // Abstract definition of a parameter value
    class Value {
    //
    public:
        // destructor.
        virtual ~Value() {}
        // Prints the parameter value to the specified stream.
        // Must be reimplemented in derived classes.
        virtual void print(std::ostream& ) const = 0;
        // Clone this value.  Useful in copy-construction.
        // Must be reimplemented in derived classes.
        virtual Value* clone() const = 0;
    };
    //*******************************************************
    // Concrete definition of a parameter value
    // for a specified type.
    template <typename T>
    class Parameter : public Value {
    //
    public:
        // \returns A read-only reference to the parameter value.
        const T& get() const { return _value; }
        // \returns A writable reference to the parameter value.
        T& set() { return _value; }
        // Prints the parameter value to the specified stream.
        virtual void print(std::ostream& ) const;
        // Clone this value.  Useful in copy-construction.
        virtual Value * clone() const;
    //
    private:
        // stored parameter value
        T _value;
    };
    //*******************************************************
    // Parameter map iterator.
    typedef std::map<std::string, Value *>::iterator iterator;
    // Constant parameter map iterator.
    typedef std::map<std::string, Value *>::const_iterator const_iterator;

    // Iterator pointing to the beginning of the set of parameters
    inline
    iterator       begin() { return _values.begin(); }
    // Iterator pointing to the beginning of the set of parameters
    inline
    const_iterator begin() const { return _values.begin(); }
    // Iterator pointing to the end of the set of parameters
    inline
    iterator       end() { return _values.end(); }
    // Iterator pointing to the end of the set of parameters
    inline
    const_iterator end() const { return _values.end(); }
//
protected:
    // data structure to map names with values
    std::map<std::string, Value*> _values;
};
// =============================================================================
// Declare this now that Parameters::print() is defined.
// By declaring this early we can use it in subsequent
// methods.  Required for gcc-4.0.2 -- 11/30/2005, BSK
inline
std::ostream& operator<<(std::ostream& os, const Parameters& p)
{
  p.print(os);
  return os;
}
// =============================================================================
// Parameters::Parameter<> class inline methods
template <typename T>
inline
void Parameters::Parameter<T>::print(std::ostream& os) const
{
  // Call helper function overloaded for basic scalar and vector types
  print_helper(os, static_cast<const T *>(&_value));
}
// -----------------------------------------------------------------------------
template <typename T>
inline
Parameters::Value* Parameters::Parameter<T>::clone() const
{
  Parameter<T>* copy = new Parameter<T>;
  copy->_value = _value;
  return copy;
}
// -----------------------------------------------------------------------------
// Parameters class inline methods
inline
void Parameters::clear() // since this is inline we must define it
{                         // before its first use (for some compilers)
  while (!_values.empty())
    {
      Parameters::iterator it = _values.begin();
      delete it->second;
      it->second = 0;
      _values.erase(it);
    }
}
// -----------------------------------------------------------------------------
inline
Parameters& Parameters::operator=(const Parameters & source)
{
  this->clear();
  *this += source;
  return *this;
}
// -----------------------------------------------------------------------------
inline
Parameters& Parameters::operator+=(const Parameters & source)
{
  Parameters::const_iterator it = source._values.begin();
  for (; it!=source._values.end(); ++it)
    {
      if (_values.find(it->first) != _values.end())
        delete _values[it->first];
      _values[it->first] = it->second->clone();
    }
  return *this;
}
// -----------------------------------------------------------------------------
inline
void Parameters::print(std::ostream& os) const
{
  Parameters::const_iterator it = _values.begin();

  os << std::endl;
  os << "Name\t Type\t Value\n"
     << "---------------------\n";
  while (it!=_values.end())
    {
      os << " " << it->first
         << "\t ";   it->second->print(os);
      os << '\n';
      //
      ++it;
    }
  os << std::flush;
}
// -----------------------------------------------------------------------------
template <typename T>
inline
bool Parameters::have_parameter(const std::string& name) const
{
  Parameters::const_iterator it = _values.find(name);

  if (it!=_values.end())
    if (dynamic_cast<const Parameters::Parameter<T>*>(it->second) != 0)
      return true;
  return false;
}
// -----------------------------------------------------------------------------
template <typename T>
inline
const T& Parameters::get(const std::string& name) const
{
  if (!this->have_parameter<T>(name))
    {
      std::ostringstream oss;
      oss << "ERROR: no";
      oss << " parameter named \""
          << name << "\" found.\n\n"
          << "Known parameters:\n"
          << *this;
      //
      ABORT_("cannot find parameter \""+name+"\"");
    }
  Parameters::const_iterator it = _values.find(name);
  return dynamic_cast<Parameter<T>*>(it->second)->get();
}
// -----------------------------------------------------------------------------
template <typename T>
inline
void Parameters::insert(const std::string & name)
{
  if (!this->have_parameter<T>(name))
    _values[name] = new Parameter<T>;
  set_attributes(name, true);
}
// -----------------------------------------------------------------------------
template <typename T>
inline
T& Parameters::set(const std::string& name)
{
  if (!this->have_parameter<T>(name))
    _values[name] = new Parameter<T>;
  set_attributes(name, false);
  return dynamic_cast<Parameter<T>*>(_values[name])->set();
}
// -----------------------------------------------------------------------------
inline
void Parameters::remove(const std::string& name)
{
  Parameters::iterator it = _values.find(name);

  if (it!=_values.end())
    {
      delete it->second;
      it->second = 0;
      _values.erase(it);
    }
}
// =============================================================================
//non-member scalar print function
template<typename P>
void print_helper(std::ostream & os, const P * param)
{
  os << (*param);
  os << std::flush;
}
// -----------------------------------------------------------------------------
template<>
inline
void print_helper(std::ostream & os, const char * param)
{
  // Specialization so that we don't print out unprintable characters
  os << static_cast<int>(*param);
  os << std::flush;
}
// -----------------------------------------------------------------------------
template<>
inline
void print_helper(std::ostream & os, const unsigned char * param)
{
  // Specialization so that we don't print out unprintable characters
  os << static_cast<int>(*param);
  os << std::flush;
}
// -----------------------------------------------------------------------------
template<typename P>
void print_helper(std::ostream & os, const std::vector<P>* param)
{
  os << std::endl;
  for (std::size_t i=0; i<param->size(); ++i)
    os << (*param)[i] << " ";
  os << std::flush;
}
// -----------------------------------------------------------------------------
template<typename P>
void print_helper(std::ostream & os, const std::vector<std::vector<P>>* param)
{
  os << std::endl;
  for (std::size_t i=0; i<param->size(); ++i)
    for (std::size_t j=0; j<(*param)[i].size(); ++j)
      os << (*param)[i][j] << " ";
  os << std::flush;
}
// =============================================================================
#endif // _MODEL_PARAMETERS_H_
// =============================================================================
