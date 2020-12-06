.. _Python render:

Rendering a scene
====================

PSDR-CUDA renders virtual scenes described using a simplified version of Mitsuba's XML `description language <https://mitsuba2.readthedocs.io/en/latest/src/getting_started/file_format.html#sec-file-format>`_.


Loading a scene
--------------------

PSDR-CUDA provides two functions for loading scenes in Python:

- ``psdr_cuda.Scene.load_file()``: load the scene from an XML file on disk;

- ``psdr_cuda.Scene.load_string()``: load the scene from a string in memory.

Here is a complete Python example on how to load a PSDR-CUDA scene from an XML file:

.. code-block:: python3

   import psdr_cuda

   scene = psdr_cuda.Scene()
   scene.load_file('scene.xml')


Forward rendering of a scene
------------------------------

Once a scene has been loaded, it can be rendered as follows:

.. code-block:: python3

   # Initialize an integrator
   integrator = psdr_cuda.DirectIntegrator()

   # Start rendering!
   sensor_id = 0
   image = integrator.renderC(scene, sensor_id)

The ``sensor_id`` variable is an integer specifying which sensor to use (as PSDR-CUDA allows a scene to contain multiple sensors with identical resolutions).

The ``image`` variable returned by the integrator's ``renderC()`` function is a RGB image stored as an Enoki CUDA array of the type ``enoki.cuda.Vector3f``. This variable can be converted to a Numpy array and saved as an OpenEXR image named ``image.exr`` using OpenCV:

.. code-block:: python3

   # Cast to a Numpy array and reshape
   out = image.numpy().reshape((scene.opts.height, scene.opts.width, 3))

   # Save to OpenEXR image using OpenCV
   import cv2
   out = cv2.cvtColor(out, cv2.COLOR_RGB2BGR) # OpenCV uses BGR
   cv2.imwrite("image.exr", out)

In the code above, ``scene.opts`` is an instance of ``psdr_cuda.RenderOption`` whose C++ definition is

.. code-block:: c++

   struct RenderOption {
       int width, height;  // Image resolution
       int spp;            // Spp for the main image/interior integral
       int sppe;           // Spp for primary edge integral
       int sppse;          // Spp for secondary edge integral
       int log_level;
   };

This variable is automatically set when loading the scene but can also be overwritten as follows.

.. code-block:: python3

   # Delay the configuration of the scene
   scene.load_file('scene.xml', auto_configure=False)

   # Overwriting the render options
   scene.opts.width     = 256
   scene.opts.height    = 256
   scene.opts.spp       = 32

   # Configure the scene
   scene.configure()

Note that a scene must be configured before being rendered. If the scene is loaded with ``auto_configure=False`` without being manually configured via ``scene.configure()``, attempting to render it will cause an assertion failure.
