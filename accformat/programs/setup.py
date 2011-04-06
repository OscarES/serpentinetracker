from distutils.core import setup, Extension
from os import environ, getcwd

uaproot     = environ.get("UAPROOT")
xercescroot = environ.get("XERCESCROOT")
pythonincs  = environ.get("PYTHONINCS")
pwd         = getcwd()

module1 = Extension('accformat',
          sources = [pwd + '/accformat.cpp'],
          include_dirs = [uaproot+'/ANTLR270',
                          uaproot,
                          xercescroot+'/src',
                          pythonincs],
          library_dirs = [uaproot+'/lib',xercescroot+'/src/.libs'],
          libraries = ['uap','antlr','xerces-c'])

setup (name = 'accformat',
       version = '1.1',
       description = 'Package for dealing with accelerator formats.',
       ext_modules = [module1])
