# lidUrb: Urban trees analyses from terrestrial laser scanning
The main functionalities are:
- **individual tree segmentation** (see ```wang_clustering``` function)
- **leaves / wood segmentation** (see ```LW_segmentation_dbscan``` and ```LW_segmentation_graph``` functions)
- **computation of leaf area density and other foliage traits** (see ```leaves_traits``` function)
- **computation of the green crown volume** (see ```green_crown_volume``` function)
- **efficient computation of points geometric features** (see ```geom_features```)
- **information and tree topology extraction from adTree QSMs** (https://github.com/tudelft3d/AdTree), (see ```adtree2qsm``` and ```qsm_topology``` functions)
- **functions to build graphs based on a point cloud** (see ```delaunay_graph```, ```knn_graph``` and ```highly_connected_graph``` functions)
- and some other functionnalities.


# Installing jakteristics

lidUrb rely on the **jakteristics** python package to compute geometric features (https://github.com/jakarto3d/jakteristics, documentation available at https://jakteristics.readthedocs.io/en/latest/) which is integrated trough the reticulate package.

You can use the following command to install it.

Check if you have python installed:
```
reticulate::py_config()
```

If python is not installed, use the ```reticulate::install_python(version)``` and select the version you want (you can use ```reticulate::install_python(list = TRUE)``` to check all available versions).

Once python is installed, you can install jakteristics using the following command:
```
reticulate::py_install("jakteristics",pip=TRUE)
```
