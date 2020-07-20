KDTree
======

.. image:: pics/nano.png

The following is a wrapper for KDTree (k-dimension tree) classes.  The main
backend that is used is `nanoflann <https://github.com/jlblancoc/nanoflann>`_,
which is a header only library that is based on the
`flann library <https://github.com/mariusmuja/flann>`_.  Nanoflann utilizes the
*curriously recurring template pattern* (CRTP) as a means of speeding up the
original flann library.

The KDTree data structure was invented in 1975 by Jon Louis Bentley [JLB]_ as a
type of binary search tree.  For a set of :math:`N` points, KDTrees have a time
complexity of :math:`\mathcal{O}(N)` for constructing the tree and an average
complexity of order :math:`\mathcal{O}(\log N)` for searching.  Worst case
scenarios for searching result in being equivalent to brute force.  As is
mentioned in the `Wikipedia <https://en.wikipedia.org/wiki/K-d_tree#High-dimensional_data>`_
article, for spaces of dimension :math:`k`, the number of points of the
data :math:`N` should be much larger than :math:`N \gg 2^k`.  Due to this problem,
many *approximate nearest neighbor* algorithms have been developed. One of these
is `FLANN <https://github.com/mariusmuja/flann>`_, which is based on the
paper by Muja and Lowe [MujaLowe]_.

Switching Backends
------------------




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
   #include <iostream>//  generate some random data
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





.. [JLB] Bentley, J. L. (1975). "Multidimensional binary search trees used for associative searching". Communications of the ACM. 18 (9): 509–517.
.. [MujaLowe] Marius Muja and David G. Lowe, "Fast Approximate Nearest Neighbors with Automatic Algorithm Configuration", in International Conference on Computer Vision Theory and Applications (VISAPP'09), 2009
