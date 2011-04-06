/*
* Program to translate an aml file and create corresponding mad8, madx, 
* xsif, and bmad lattice files.
*/

#include "Python.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "UAP/UAPUtilities.hpp"
#include "Translate/MAD8Parser.hpp"
#include "Translate/MADXParser.hpp"
#include "Translate/BmadParser.hpp"
#include "Translate/XSIFParser.hpp"
#include "Translate/AMLParser.hpp"
#include "Translate/SADParser.hpp"
#include "AML/AMLReader.hpp"
#include "AML/AMLLatticeExpander.hpp"

using namespace std;

static PyObject *
convaml2xsif (PyObject *self, PyObject *args) { 

  UAPNode* UAPModel = NULL;
  const char *file_name;

  if (!PyArg_ParseTuple(args, "s", &file_name))
      return NULL;

  // Read in the AML file.
  AMLReader reader;
  UAPModel = reader.AMLFileToAMLRep (file_name);
  if (!UAPModel) UAPModel = reader.AMLFileToAMLRep(((string)file_name).append(".aml"));
  if (!UAPModel)
  {
    PyErr_Print();
    PyErr_SetString(PyExc_IOError, "File does not exist.");
    return NULL;
  }

  // Convert to MAD format, etc. and output a file for each.
  // First strip off ".aml" if the input file has this suffix
  string file_stub = file_name;
  size_t ix_aml = file_stub.find(".aml");
  if (ix_aml != string::npos) file_stub.erase(ix_aml);

  TranslateCore* writer = new XSIFParser();
  string out_file = file_stub + "_out.xsif";
  writer->AMLRepToXFile(UAPModel, out_file);
  delete writer;

  return Py_BuildValue("i", 1);
}

static PyObject *
convaml2madx (PyObject *self, PyObject *args) {

  UAPNode* UAPModel = NULL;
  const char *file_name;

  if (!PyArg_ParseTuple(args, "s", &file_name))
      return NULL;

  // Read in the AML file.
  AMLReader reader;
  UAPModel = reader.AMLFileToAMLRep (file_name);
  if (!UAPModel) UAPModel = reader.AMLFileToAMLRep(((string)file_name).append(".aml"));
  if (!UAPModel)
  {
    PyErr_Print();
    PyErr_SetString(PyExc_IOError, "File does not exist.");
    return NULL;
  }

  // Convert to MAD format, etc. and output a file for each.
  // First strip off ".aml" if the input file has this suffix
  string file_stub = file_name;
  size_t ix_aml = file_stub.find(".aml");
  if (ix_aml != string::npos) file_stub.erase(ix_aml);

  TranslateCore* writer = new MADXParser();
  string out_file = file_stub + ".out.madx";
  writer->AMLRepToXFile(UAPModel, out_file);
  delete writer;

  return Py_BuildValue("i", 1);
}

static PyObject *
convaml2mad8 (PyObject *self, PyObject *args) {

  UAPNode* UAPModel = NULL;
  const char *file_name;

  if (!PyArg_ParseTuple(args, "s", &file_name))
      return NULL;

  // Read in the AML file.
  AMLReader reader;
  UAPModel = reader.AMLFileToAMLRep (file_name);
  if (!UAPModel) UAPModel = reader.AMLFileToAMLRep(((string)file_name).append(".aml"));
  if (!UAPModel)
  {
    PyErr_Print();
    PyErr_SetString(PyExc_IOError, "File does not exist.");
    return NULL;
  }

  // Convert to MAD format, etc. and output a file for each.
  // First strip off ".aml" if the input file has this suffix

  string file_stub = file_name;
  size_t ix_aml = file_stub.find(".aml");
  if (ix_aml != string::npos) file_stub.erase(ix_aml);

  TranslateCore* writer = new MAD8Parser();
  string out_file = file_stub + ".out.mad8";
  writer->AMLRepToXFile(UAPModel, out_file);
  delete writer;

  return Py_BuildValue("i", 1);
}

static PyObject *
convaml2bmad (PyObject *self, PyObject *args) {

  UAPNode* UAPModel = NULL;
  const char *file_name;

  if (!PyArg_ParseTuple(args, "s", &file_name))
      return NULL;

  // Read in the AML file.
  AMLReader reader;
  UAPModel = reader.AMLFileToAMLRep (file_name);
  if (!UAPModel) UAPModel = reader.AMLFileToAMLRep(((string)file_name).append(".aml"));
  if (!UAPModel)
  {
    PyErr_Print();
    PyErr_SetString(PyExc_IOError, "File does not exist.");
    return NULL;
  }

  // Convert to MAD format, etc. and output a file for each.
  // First strip off ".aml" if the input file has this suffix
  string file_stub = file_name;
  size_t ix_aml = file_stub.find(".aml");
  if (ix_aml != string::npos) file_stub.erase(ix_aml);

  TranslateCore* writer = new BmadParser();
  string out_file = file_stub + ".out.bmad";
  writer->AMLRepToXFile(UAPModel, out_file);
  delete writer;

  return Py_BuildValue("i", 1);
}

static PyObject *
amlreader (PyObject *self, PyObject *args){
  UAPNode *UAPModel = NULL, *MLattice = NULL;
  const char *file_name, *outstring;
  string accstringrep;

  if (!PyArg_ParseTuple(args, "s", &file_name)) return NULL;

  // Read in the AML file.
  AMLReader reader;
  UAPModel = reader.AMLFileToAMLRep (file_name);
  if (!UAPModel) UAPModel = reader.AMLFileToAMLRep(((string)file_name).append(".aml"));
  if (!UAPModel) {
    PyErr_Print();
    PyErr_SetString(PyExc_IOError, "File does not exist.");
    return NULL;
  }

  AMLLatticeExpander LE;
  LE.AMLExpandLattice(UAPModel);

  MLattice = (UAPModel->getChildByName("expanded_lattice"))->getChildByName("machine");
  if (!MLattice) return NULL;

  accstringrep = MLattice->toXMLTree();
  outstring = accstringrep.c_str();
  return Py_BuildValue("s", outstring);
}

static PyObject *
latticereader (PyObject *self, PyObject *args){
  UAPNode *UAPModel = NULL, *MLattice = NULL;
  const char *fnamecnstchar, *outstring;
  string file_name,accstringrep;

  if (!PyArg_ParseTuple(args, "s", &fnamecnstchar)) return NULL;
  file_name = fnamecnstchar;

  TranslateCore* reader;
  if      (file_name.find(".madx") != string::npos) reader = new MADXParser();
  else if (file_name.find(".bmad") != string::npos) reader = new BmadParser();
  else if (file_name.find(".mad8") != string::npos) reader = new MAD8Parser();
  else if (file_name.find(".xsif") != string::npos) reader = new XSIFParser();
  else if (file_name.find(".aml")  != string::npos) reader = new AMLParser();
  else if (file_name.find(".sad")  != string::npos) reader = new SADParser();
  else {
    PyErr_Print();
    PyErr_SetString(PyExc_IOError, "Error: Unknown file extension for input file");
    return NULL;
  }

  UAPModel = reader->XFileToAMLRep (file_name);
  if (!UAPModel) {
    PyErr_Print();
    PyErr_SetString(PyExc_IOError, "File does not exist.");
    return NULL;
  }

  AMLLatticeExpander LE;
  LE.AMLExpandLattice(UAPModel);

  MLattice = (UAPModel->getChildByName("expanded_lattice"))->getChildByName("machine");
  if (!MLattice) return NULL;

  accstringrep = MLattice->toXMLTree();
  outstring = accstringrep.c_str();
  return Py_BuildValue("s", outstring);
}

static PyMethodDef AccformatMethods[] = {
  {"convaml2xsif", convaml2xsif, METH_VARARGS, "Convert AML file to XSIF representation."},
  {"convaml2madx", convaml2madx, METH_VARARGS, "Convert AML file to MADX representation."},
  {"convaml2mad8", convaml2mad8, METH_VARARGS, "Convert AML file to MAD8 representation."},
  {"convaml2bmad", convaml2bmad, METH_VARARGS, "Convert AML file to BMAD representation."},
  {"amlreader", amlreader, METH_VARARGS, "Read an AML file, and convert it to a heirarchical string representation."},
  {"latticereader", latticereader, METH_VARARGS, "Read a lattice file, and convert it to a heirarchical string representation."},
  {NULL, NULL, 0, NULL}    /* Sentinel */
};

PyMODINIT_FUNC
initaccformat(void)
{
  (void) Py_InitModule("accformat", AccformatMethods);
}

