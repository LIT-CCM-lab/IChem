# Python API - Install

 
This page explains how to use the **Python bindings** (`ichem_ifp`).


There are two common ways to get the module:

1. **Recommended:** download a **prebuilt IChem bundle** (Ubuntu/Debian) that includes the Python extension
2. **Developer:** build the extension from source with CMake

---

## Option 1: Install from a prebuilt bundle (recommended)

IChem bundles ship both:
- the CLI binary: `bin/IChem`
- the Python extension: `python/ichem_ifp.cpython-312-*.so`

### 1) Extract the archive

```bash
tar -xzf ichem_ubuntu_x86_64_v5.3.8.tar.gz
cd ichem_ubuntu_x86_64_v5.3.8
```

### 2) Import the module (no installation required)
The module is located in python/. Add it to PYTHONPATH:

PYTHONPATH="$(pwd)/python" python3.12 -c "import ichem_ifp; print(ichem_ifp.__file__)"

### 3) Run your script
export PYTHONPATH="$(pwd)/python:$PYTHONPATH"
python3.12 test_ifp_api.py

Note: The extension filename contains cpython-312, so you must use Python 3.12 to import it.

## Option 2 —
- Install with pip (wheel)

( need to ship many versions, now only 3.12 and over)
provide a .whl compatible with Python 3.12:
IChem ships prebuilt Python wheels. Install the wheel matching your Python version and platform.

1) Download according to your target version
2) pip install ./ichem_ifp-*.whl


## Option 3: Build from source (developer)

For full instructions, see the root `README.md`.
