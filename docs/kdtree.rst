KDTree
======

.. image:: pics/nano.png

The following is a wrapper for KDTree (k-dimension tree) classes.  The main
backend that is used is `nanoflann <https://github.com/jlblancoc/nanoflann>`_,
which is a header only library that is based on the
`flann library <https://github.com/mariusmuja/flann>`_.  Nanoflann utilizes the
*curriously recurring template pattern* (CRTP) as a means of speeding up the
original flann library.

Creating a KDTree
-----------------

A KDTree can be created by simply calling the default constructor with
an appropriate template parameter.  The following shows an example for
a random two-dimensional data set.

.. code-block:: c
   :linenos:

   //  include the kdtree header
   #include "kdtree.h"
   #include <random>
   #include <iostream>
   int main
   {
      //  generate some random data
      const int range_from  = 0.0;
      const int range_to    = 1.0;
      std::random_device rand_dev;
      std::mt19937 generator(rand_dev());
      std::uniform_real_distribution<double> distr(range_from, range_to);

      size_t N = 10000;
      std::vector<std::vector<double>> data(N);
      for (auto i = 0; i < N; i++) {
        data[i][0] = distr(generator);
        data[i][1] = distr(generator);
      }

      ET::KDTree<double> kdt = new ET::KDTree<double>(data);

      //  Query k neighbors for each point in data
      size_t k = 5
      kdt->queryNeighbors(k);

      delete kdt;
      return 0;
   }

Here is a similar example using the python bindings.

.. code-block:: python
   :linenos:

   # import the KDTree class
   from etraj import KDTree
   import numpy as np

   low = 0.0
   high = 1.0
   N = 10000
   data = np.random.uniform(low,high,(N,2))

   # instantiate the tree with the data
   kdt = KDTree(data)

   # Query k neighbors for each point in data
   k = 5
   kdt.query_neighbors(k)






KDTree Source Code
------------------

.. doxygenclass:: ET::KDTree
   :project: etraj
   :members:
