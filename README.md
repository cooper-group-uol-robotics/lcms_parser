# Parser for UPLC-MS raw data in A.I.C. Group
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11174536.svg)](https://doi.org/10.5281/zenodo.11174536)

## Requirements and installation

Requires the official Python SDK from Waters. On runtime, a valid license will also be needed.

### MassLynxSDK

The [MassLynxSDK](https://microapps.on-demand.waters.com/home/downloads/masslynx-sdk) needs to be requested from Waters directly and the EULA will need to be accepted. It can be installed within the desired environment directly from the wheel file provided with the package - it can be found in one of the zipped directories:

```
pip install masslynxsdk-5.0.0-py3-none-any.whl
```

### MassLynxSDK license

The license key (string found in the `license.key` file) needs to be provided when interfacing with the API:

```python
from lcms_parser import WatersRawFile

raw_file = WatersRawFile(path="test.raw", license_key="XXXXXX")
```

Alternatively, the `license_key` parameter can be omitted as long as a valid `license.key` file is present in the current working directory.

## Notes to developers and contributors

There are linters and formatters in use for this project. Prior to contributing code, please make sure that your development environment is set up. Typically, an editable version would be installed for development and `pre-commit` would ensure that the code conforms to the standards before it is commmitted to GitHub:

```
python -m venv venv
venv\Scripts\activate
pip install -e .[dev] --find-links wheels
pre-commit install
```

Ensure that the [`masslynxsdk`](https://microapps.on-demand.waters.com/home/downloads/masslynx-sdk) from Waters is either installed on the system or that the wheel is present in `wheels/`.