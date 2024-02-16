# pyNorm: installing / using `pyNorm`

There are several ways to make use of `pyNorm`. Currently the principal approaches are:

## **Editable pip installation:** 

This installation approach allows edits to the code / `git pull` updates to be directly accessible. To enable this installation, invoke the following in the terminal from within the `pyNorm` directory:

```
pip install -e .
```

## **Install from GitHub via `pip`:**

```
pip install git+https://github.com/jchowk/pyNorm.git
```

## **Include `pyNorm` in your `$PYTHONPATH`:**

This approach also makes changes to the base code immediately accessible. Add the full path to the `pyNorm` code to your `$PYTHONPATH` variable by invoking something like (from the `bash` terminal or within your shell configuration file):

```
export PYTHONPATH="$PYTHONPATH:/path/to/pyNorm/pyNorm/"
```

Note the path has to point to the subdirectory `pyNorm/pyNorm/`. 
