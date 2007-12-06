// Patrick Hung
// California Inst. of Tech.


#include <portinfo>
#include <Python.h>
// hack for Python 2.3
#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE return Py_INCREF(Py_None), Py_None
#endif

#include "phonexterns.h"

#include "journal/debug.h"

#include <vector>
#include <iostream>

#ifdef USE_NUMERIC
#include "Numeric/arrayobject.h"
#endif
#ifdef USE_NUMPY
#include "numpy/arrayobject.h"
#endif

#define JINFO info << journal::at(__HERE__)

//extern "C" void Py_DebugTrap(void);

static PyObject *_pyphonError;

PyDoc_STRVAR(
    module_doc,
    "_pyphon \n"
    "\n"
    "    -- Python extension module to  Finite Element Method on Pyre\n"
    "\n"
    "\n"
    );


PyDoc_STRVAR(hello_doc, "hello world for the module");
static PyObject * _pyphon_hello(PyObject *self, PyObject *args) {
    return PyString_FromString("Hello, c'est moi you're looking for ?");
}

PyDoc_STRVAR(phon_doc, "Run phon (the executable)");
static PyObject * _pyphon_phon(PyObject *self, PyObject *args) {
    phonfunc(); 
    return PyString_FromString("XXX");
    //Py_RETURN_NONE;
}

PyDoc_STRVAR(hpsort_doc, "Wrapper for hpsort.");
static PyObject * _pyphon_hpsort(PyObject *self, PyObject *args) {
    PyObject *in;
    int n;

    if (! PyArg_ParseTuple(args, "O:hpsort", &in)) return 0;

    Py_INCREF(in);

    if (!PySequence_Check(in)) {
        return 0;
    }
    n = PySequence_Length(in);
    PyObject * py_floats = PySequence_Fast(in, "Expecting a sequence of floats.");

    std::vector<double> tobesorted(n);
    std::vector<int> index(n);

    for (int i=0; i<n; ++i) {
        tobesorted[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(py_floats, i));
    }
    Py_DECREF(in);

    hpsort(&n, &tobesorted[0], &index[0]);

    PyArrayObject *n1, *n2;
    int dims[] = {n};
    int nd = 1;
    n1 = (PyArrayObject*) PyArray_FromDims(nd, dims, PyArray_DOUBLE);
    n2 = (PyArrayObject*) PyArray_FromDims(nd, dims, PyArray_INT);
    const int s0 = n1->strides[0];
    const int s1 = n2->strides[0];
    for (int i=0; i<n; ++i) {
        double *xx = (double*) (n1->data + i*s0);
        int *yy = (int*) (n2->data + i*s1);
        xx[0] = tobesorted[i];
        yy[0] = index[i];
    }

    PyObject *ret;
    ret = PyTuple_New(2);
    PyTuple_SetItem(ret, 0, PyArray_Return(n1));
    PyTuple_SetItem(ret, 1, PyArray_Return(n2));
    return ret;
}


// List of functions defined in this module
struct PyMethodDef _pyphon_methods[] = {

    // the core routines
    { "hello",   _pyphon_hello,   1, hello_doc },
    { "phon",   _pyphon_phon,   1, phon_doc },
    { "hpsort",   _pyphon_hpsort,   1, hpsort_doc },

    {  NULL, NULL } // sentinel

};


PyMODINIT_FUNC init_pyphon(void)
{
    PyObject * m;

    // Create the module and add the functions
    m = Py_InitModule3( "_pyphon", _pyphon_methods, module_doc);

    _pyphonError = PyErr_NewException("_pyphon.error", NULL, NULL);
    Py_INCREF(_pyphonError);
    PyModule_AddObject(m, "error", _pyphonError);

    import_array();

    // check for errors
    if (PyErr_Occurred()) {
        Py_FatalError("can't initialize module _pyphon");
    }

    return;
}

#undef JINFO 

// end of file
