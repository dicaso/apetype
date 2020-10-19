.. apetype documentation master file, created by
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to APEtype's documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   configurations
   tasks
   multiprocessing
   reports
   examples

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

  
Introduction
============

"To type or not to type?", that is the age-old question in
programming. To avoid a steep learning curve and make explorative
programming fun, python has opted for dynamic typing. Recently type
hints have been introduced, to accomodate a need for stronger typing
context. This gap in python typing, had been partly met in the past by
packages such as ``luigi``, developed at *Spotify*.

I used ``luigi`` as the base for my pipelines and bioinformatics
workflows in my tool ``genairics`` before being aware of the type
hints. As I learned about the type hints, I quickly realised the
``luigi.Parameter`` could be more elegantly replaced by type hints.

``apetype`` takes the concept of ``luigi.Parameter`` and function
injections by the TypeScript ``Angular`` platform, to create a new
library for creating workflows and pipelines, that can easily be
started as function calls or command line executables out of the
box. The workflows can also divide their tasks across the multiple
processors in a system, and in the future, with few modifications of
the current code base, they should also be deployable on the cloud.
  
Quick start
===========
.. automodule:: apetype
   :members:

