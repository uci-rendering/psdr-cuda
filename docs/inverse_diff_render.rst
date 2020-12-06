.. _Inverse diff render:

Differentiable rendering
=========================

How to use PSDR-CUDA for forward rendering has been shown in the :ref:`previous section <Python render>`.
In what follows, we focus on performing *differentiable rendering*.

Similar to the forward-rendering case, after being loaded and configured, a scene can be rendered in a differentiable fashion by calling an integrator's ``renderD()`` method. The return value of this method is an Enoki CUDA array of the type ``enoki.cuda_autodiff.Vector3f``.


Generating derivative image
---------------------------------------

The following example generates a derivative image ``image_deriv`` with respect to the rotation angle ``P`` of the first mesh about the Z-axis.

.. code-block:: python3

   import enoki as ek
   from enoki.cuda_autodiff import Float32 as FloatD, Vector3f as Vector3fD, Matrix4f as Matrix4fD
   import psdr_cuda

   # Load the scene without configuring it
   scene = psdr_cuda.Scene()
   scene.load_file('scene.xml', auto_configure=False)

   # Initialize a parameter P for differentiation
   P = FloatD(0.)
   ek.set_requires_gradient(P)

   # Construct a rotation matrix
   mat = Matrix4fD.rotate(Vector3fD(0., 0., 1.), P)

   # Left-multiply this matrix to the to-world transformation of the first mesh of the scene
   scene.param_map["Mesh[0]"].append_transform(mat)

   # Configure the scene
   scene.configure()

   # Start rendering!
   image = psdr_cuda.DirectIntegrator().renderD(scene, sensor_id=0)

   # Compute derivative image with respect to P using forward-mode autodiff
   ek.forward(P)
   image_deriv = ek.gradient(image)


Differentiating image loss
---------------------------------------

What follows is another example that computes the gradient of the image RMSE ``loss``, given a target image ``target_image``, with respect to all vertex positions of ``Mesh[0]``. The resulting ``grad`` can then be used in gradient-based optimization pipelines to minimize ``loss``.

.. code-block:: python3

   import enoki as ek
   from enoki.cuda_autodiff import Float32 as FloatD, Vector3f as Vector3fD, Matrix4f as Matrix4fD
   import psdr_cuda

   # Load the scene without configuring it
   scene = psdr_cuda.Scene()
   scene.load_file('scene.xml', auto_configure=False)

   # Compute gradient with respect to mesh vertex positions
   ek.set_requires_gradient(scene.param_map["Mesh[0]"].vertex_positions)

   # Configure the scene
   scene.configure()

   # Start rendering!
   image = psdr_cuda.DirectIntegrator().renderD(scene, sensor_id=0)

   # Compute the RMSE image loss
   loss = ek.sqrt(ek.hmean(ek.squared_norm(target_image - image)))

   # Reverse-mode autodiff
   ek.backward(loss)

   # Obtain the gradient of the loss
   grad = ek.gradient(scene.param_map["Mesh[0]"].vertex_positions)
