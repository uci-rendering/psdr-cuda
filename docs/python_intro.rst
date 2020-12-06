.. _Python intro:

Introduction
====================

Although the core functionalities of PSDR-CUDA are written in C++17 and CUDA, to enable easy integration into Python-based inference and learning pipelines, running PSDR-CUDA requires using its Python bindings.

This section focuses on the use of the bindings for *forward* rendering.
We will show how to perform differentiable and inverse rendering in the :ref:`next section <Inverse intro>`.


Importing the renderer
------------------------------

To import PSDR-CUDA into Python, you only need to import its central module named ``psdr_cuda``:

.. code-block:: python3

   import psdr_cuda


Basic arithmetic types
------------------------------

PSDR-CUDA relies heavily on the `Enoki <https://github.com/mitsuba-renderer/enoki/>`_ library for elementary arithmetic types and mathematical operations. Commonly used data types include ``Float32``, ``Vector3f``, and ``Matrix4f`` from both ``enoki.cuda``, which can be imported as follows:

.. code-block:: python3

   from enoki.cuda import Float32 as FloatC, Vector3f as Vector3fC, Matrix4f as Matrix4fC

For more details on Enoki's data types and their Python bindings, please refer to its `official documentation <http://enoki.readthedocs.org/en/master>`_.
