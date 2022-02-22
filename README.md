# lidUrb
Urban trees analyses from terrestrial laser scanning

# Installing jakteristics

lidUrb partially rely on the **jakteristics** python package to compute geometric features (https://github.com/jakarto3d/jakteristics, documentation available at https://jakteristics.readthedocs.io/en/latest/) which is integrated trough the reticulate package.

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
