.. _Inverse intro:

Introduction
====================

We have showed how to perform forward rendering using PSDR-CUDA in the :ref:`previous section <Python intro>`.
This section focuses on the use of PSDR-CUDA's Python bindings for *differentiable* and *inverse* rendering applications.


Automatic-differentiation types
-----------------------------------

A key ingredient to differentiable rendering is automatic differentiation. PSDR-CUDA relies on Enoki's to perform transparent reverse-mode automatic differentiation. To this end, commonly used data types include ``Float32``, ``Vector3f``, and ``Matrix4f`` from ``enoki.cuda_autodiff``, which can be imported as follows:

.. code-block:: python3

   from enoki.cuda_autodiff import Float32 as FloatD, Vector3f as Vector3fD, Matrix4f as Matrix4fD

For more details on how Enoki's automatic differentiation works, please refer to its `official documentation <http://enoki.readthedocs.org/en/master/autodiff.html>`_.
